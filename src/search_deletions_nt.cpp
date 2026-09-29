/*
 * This File is part of Pindel; a program to locate genomic variation.
 * https://trac.nbic.nl/pindel/
 *
 *   Copyright (C) 2011 Kai Ye
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// Pindel header files
#include "logstream.h"
#include "reporter.h"
#include "control_state.h"
#include "logdef.h"

#include <string>
#include <vector>

unsigned editDistance(const std::string& a, const std::string& b)
{
   const unsigned n = a.size();
   const unsigned m = b.size();
   std::vector<unsigned> prev(m + 1), cur(m + 1);
   for (unsigned j = 0; j <= m; ++j) {
      prev[j] = j;
   }
   for (unsigned i = 1; i <= n; ++i) {
      cur[0] = i;
      for (unsigned j = 1; j <= m; ++j) {
         const unsigned cost = (a[i - 1] == b[j - 1]) ? 0u : 1u;
         unsigned best = prev[j] + 1;
         if (cur[j - 1] + 1 < best) {
            best = cur[j - 1] + 1;
         }
         if (prev[j - 1] + cost < best) {
            best = prev[j - 1] + cost;
         }
         cur[j] = best;
      }
      prev.swap(cur);
   }
   return prev[m];
}

/* A long anchor gap that still matches the read, apart from a few inserted
   bases and scattered mismatches, is a short insertion. A mismatch at either
   end must not keep the whole span as one replacement. */
bool gapIsShortInsertion(const std::string& chromosome, const SPLIT_READ& read)
{
   if (read.IndelSize == 0 || read.NT_size <= read.IndelSize || read.NT_str.size() != read.NT_size) {
      return false;
   }
   const unsigned extra = read.NT_size - read.IndelSize;
   if (extra > 5) {
      return false;
   }
   const unsigned start = (unsigned) read.Left + (unsigned) read.BP + 1;
   if (start + read.IndelSize > chromosome.size()) {
      return false;
   }
   const std::string refGap = chromosome.substr(start, read.IndelSize);
   const unsigned allowed = extra + 2 + refGap.size() / 5;
   return editDistance(refGap, read.NT_str) <= allowed;
}

/* Place `extra` read bases into the reference gap. The best match, and the
   leftmost one on a tie, is where the short insertion sits. */
bool locateInsertion(const std::string& refGap, const std::string& nt,
                     unsigned extra, unsigned& refOff, unsigned& ntOff)
{
   if (extra == 0 || nt.size() != refGap.size() + extra) {
      return false;
   }
   bool found = false;
   unsigned bestScore = 0;
   unsigned bestAt = 0;
   for (unsigned insAt = 0; insAt + extra <= nt.size() && insAt <= refGap.size(); ++insAt) {
      unsigned score = 0;
      for (unsigned i = 0; i < refGap.size(); ++i) {
         const unsigned j = (i < insAt) ? i : i + extra;
         if (j < nt.size() && refGap[i] == nt[j]) {
            ++score;
         }
      }
      if (!found || score > bestScore) {
         found = true;
         bestScore = score;
         bestAt = insAt;
      }
   }
   if (!found || bestScore * 5 < refGap.size() * 4) {
      return false;
   }
   refOff = bestAt;
   ntOff = bestAt;
   return true;
}

/* Point the anchors at the inserted base so the short-insertion search can
   aggregate every read there, instead of keeping each mismatched span. */
bool retargetShortInsertion(SPLIT_READ& read, const std::string& chromosome)
{
   const unsigned extra = read.NT_size - read.IndelSize;
   const unsigned gapStart = (unsigned) read.Left + (unsigned) read.BP + 1;
   if (gapStart + read.IndelSize > chromosome.size() || gapStart == 0) {
      return false;
   }
   const std::string refGap = chromosome.substr(gapStart, read.IndelSize);
   unsigned refOff = 0;
   unsigned ntOff = 0;
   if (!locateInsertion(refGap, read.NT_str, extra, refOff, ntOff)) {
      return false;
   }
   const int leftLen = (int) read.BP + 1 + (int) ntOff;
   const int rightLen = read.getReadLength() - leftLen - (int) extra;
   if (leftLen < 1 || rightLen < 1 || leftLen > Max_short || rightLen > Max_short) {
      return false;
   }
   const unsigned leftAbs = gapStart + refOff - 1;
   const unsigned rightAbs = leftAbs + 1;
   if (read.MatchedD == Plus) {
      const UniquePoint oldClose = read.UP_Close[read.UP_Close.size() - 1];
      const UniquePoint oldFar = read.UP_Far[read.UP_Far.size() - 1];
      read.UP_Close.clear();
      read.UP_Far.clear();
      read.UP_Close.push_back(UniquePoint(oldClose.chromosome_p, (short) leftLen, leftAbs,
                                           oldClose.Direction, oldClose.Strand, 0));
      read.UP_Far.push_back(UniquePoint(oldFar.chromosome_p, (short) rightLen, rightAbs,
                                         oldFar.Direction, oldFar.Strand, 0));
   } else if (read.MatchedD == Minus) {
      const UniquePoint oldClose = read.UP_Close[read.UP_Close.size() - 1];
      const UniquePoint oldFar = read.UP_Far[read.UP_Far.size() - 1];
      read.UP_Close.clear();
      read.UP_Far.clear();
      read.UP_Far.push_back(UniquePoint(oldFar.chromosome_p, (short) leftLen, leftAbs,
                                         oldFar.Direction, oldFar.Strand, 0));
      read.UP_Close.push_back(UniquePoint(oldClose.chromosome_p, (short) rightLen, rightAbs,
                                           oldClose.Direction, oldClose.Strand, 0));
   } else {
      return false;
   }
   return true;
}

int searchIndels(ControlState& currentState, unsigned NumBoxes, const SearchWindow& window )
{

   static int Count_DI = 0;
   static int Count_DI_Plus = 0;
   static int Count_DI_Minus = 0;

   unsigned CloseIndex, FarIndex;

   std::vector<unsigned> DI[NumBoxes];
   unsigned TempBoxIndex;
   LOG_INFO(*logStream << "Searching deletion-insertions ... " << std::endl);
   const std::string& chromosomeSeq = window.getChromosome()->getSeq();

   //UserDefinedSettings *userSettings = UserDefinedSettings::Instance();

   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++) {
      SPLIT_READ& currentRead = currentState.Reads_SR[ReadIndex];
      if (currentRead.Used
            || currentRead.UP_Far.empty() || currentRead.FragName != currentRead.FarFragName) {
         continue;
      }

      CloseIndex = currentRead.UP_Close.size() - 1;
      FarIndex = currentRead.UP_Far.size() - 1;
      if (currentRead.UP_Far[FarIndex].Mismatches + currentRead.UP_Close[CloseIndex].Mismatches
            > (short) (1 + userSettings->Seq_Error_Rate * (currentRead.UP_Far[FarIndex].LengthStr + currentRead.UP_Close[CloseIndex].LengthStr))) {
         continue;
      }

      if (currentRead.MatchedD == Plus) {
         if (currentRead.UP_Far[FarIndex].Direction == Minus) {
            if (currentRead.UP_Far[FarIndex].LengthStr + currentRead.UP_Close[CloseIndex].LengthStr < currentRead.getReadLength() &&
                  currentRead.UP_Far[FarIndex].LengthStr + currentRead.UP_Close[CloseIndex].LengthStr >= userSettings->Min_Num_Matched_Bases &&
                  currentRead.UP_Far[FarIndex].AbsLoc > currentRead.UP_Close[CloseIndex].AbsLoc + 1) {
               currentRead.Left = currentRead.UP_Close[CloseIndex].AbsLoc - currentRead.UP_Close[CloseIndex].LengthStr + 1;
               currentRead.Right = currentRead.UP_Far[FarIndex].AbsLoc + currentRead.UP_Far[FarIndex].LengthStr - 1;
               currentRead.BP = currentRead.UP_Close[CloseIndex].LengthStr - 1;
               currentRead.NT_size = currentRead.getReadLength() - currentRead.UP_Far[FarIndex].LengthStr - currentRead.UP_Close[CloseIndex].LengthStr;

               currentRead.NT_str = ReverseComplement( currentRead. getUnmatchedSeq()). substr( currentRead.BP + 1, currentRead.NT_size);
               currentRead.IndelSize = (currentRead.Right - currentRead.Left) + currentRead.NT_size - currentRead.getReadLengthMinus();

               currentRead.BPLeft = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
               currentRead.BPRight = currentRead.UP_Far[FarIndex].AbsLoc - g_SpacerBeforeAfter;
               if (gapIsShortInsertion(chromosomeSeq, currentRead)) {
                  retargetShortInsertion(currentRead, chromosomeSeq);
                  continue;
               }

               if (1) {

                  if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                     saveReadForNextCycle(currentRead, currentState.FutureReads_SR);
                  } else {
                     if (readInSpecifiedRegion( currentRead, userSettings->getRegion() ) ) {
                        TempBoxIndex = (int) (currentRead. BPLeft) / BoxSize;
                        if (TempBoxIndex < NumBoxes) {

                           DI[TempBoxIndex]. push_back(ReadIndex);
                           currentRead.Used = true;
                           Count_DI++;
                           Count_DI_Plus++;
                        }
                     }
                  }
               }
            }
         }
      } else if (currentRead.MatchedD == Minus) {
         if (currentRead.UP_Far[FarIndex].Direction == Plus) {
            if (currentRead.UP_Close[CloseIndex].LengthStr + currentRead.UP_Far[FarIndex].LengthStr < currentRead.getReadLength() &&
                  currentRead.UP_Close[CloseIndex].LengthStr + currentRead.UP_Far[FarIndex].LengthStr >= userSettings->Min_Num_Matched_Bases &&
                  currentRead.UP_Close[CloseIndex].AbsLoc > currentRead.UP_Far[FarIndex].AbsLoc + 1) {

               currentRead.Left = currentRead.UP_Far[FarIndex].AbsLoc - currentRead.UP_Far[FarIndex].LengthStr + 1;
               currentRead.Right = currentRead.UP_Close[CloseIndex].AbsLoc + currentRead.UP_Close[CloseIndex].LengthStr - 1;
               currentRead.BP = currentRead.UP_Far[FarIndex].LengthStr - 1;
               currentRead.NT_size = currentRead.getReadLength() - currentRead.UP_Close[CloseIndex].LengthStr - currentRead.UP_Far[FarIndex].LengthStr;
               currentRead.NT_str = currentRead.getUnmatchedSeq(). substr( currentRead.BP + 1, currentRead.NT_size);

               currentRead.IndelSize = (currentRead.Right - currentRead.Left) - currentRead.getReadLengthMinus() + currentRead.NT_size;
               currentRead.BPLeft = currentRead.UP_Far[FarIndex].AbsLoc - g_SpacerBeforeAfter;
               currentRead.BPRight = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
               if (gapIsShortInsertion(chromosomeSeq, currentRead)) {
                  retargetShortInsertion(currentRead, chromosomeSeq);
                  continue;
               }
               {
                  if ( 1 ) {

                     if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                        saveReadForNextCycle( currentRead, currentState.FutureReads_SR);
                     } else {
                        if (readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                           TempBoxIndex = (int) (currentRead. BPLeft) / BoxSize;
                           if (TempBoxIndex < NumBoxes) {

                              DI[TempBoxIndex]. push_back(ReadIndex);
                              currentRead.Used = true;
                              Count_DI++;
                              Count_DI_Minus++;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   LOG_INFO(*logStream << "Total: " << Count_DI << "\t+" << Count_DI_Plus << "\t-"
            << Count_DI_Minus << std::endl);
   std::ofstream DeletionOutf( userSettings->getDOutputFilename().c_str(), std::ios::app);
   std::ofstream inversionsOutf( userSettings->getINVOutputFilename().c_str(), std::ios::app);
   SortOutputDI(currentState, NumBoxes, window.getChromosome()->getSeq(), currentState.Reads_SR, DI, DeletionOutf, inversionsOutf);
   DeletionOutf.close();
   for (unsigned int i = 0; i < NumBoxes; i++) {
      DI[i].clear();
   }

   return EXIT_SUCCESS;
}
