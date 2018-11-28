#!./lib/bats/bin/bats

load 'lib/bats-support/load'
load 'lib/bats-assert/load'

load test-helper

#PINDEL=~/tmp/pindel/pindel
PINDEL=../../pindel

#PINDEL2VCF=~/tmp/pindel/pindel2vcf
PINDEL2VCF=../../pindel2vcf

DEMODIR=../../demo
DEMOREF=$DEMODIR/simulated_reference.fa
REALREF=hs_ref_chr5.fa

# --------------------------------------------------------------------------------------
#
#                             p i n d e l ,   d e m o
#
@test "pindel: demo BAMs ➔ D, TD, INV, SI" {
  run $PINDEL -i simulated_config.txt -f $DEMOREF -o result/no-strand-counts -c ALL
  assert_success
}

@test "pindel: demo D" {
  run cmp result/no-strand-counts_D expected/no-strand-counts_D
  assert_success
}

@test "pindel: demo TD" {
  run cmp result/no-strand-counts_TD expected/no-strand-counts_TD
  assert_success
}

@test "pindel: demo INV" {
  run cmp result/no-strand-counts_INV expected/no-strand-counts_INV
  assert_success
}

@test "pindel: demo SI" {
  run cmp result/no-strand-counts_SI expected/no-strand-counts_SI
  assert_success
}


# --------------------------------------------------------------------------------------
#
#                             p i n d e l ,   r e a l
#
@test "pindel: real BAMs ➔ D, INV" {
  run $PINDEL -i real_config.txt -f $REALREF -o result/no-strand-counts-real -c ALL
  assert_success
}

@test "pindel: real D" {
  run cmp result/no-strand-counts-real_D expected/no-strand-counts-real_D
  assert_success
}

@test "pindel: real INV" {
  run cmp result/no-strand-counts-real_INV expected/no-strand-counts-real_INV
  assert_success
}


# --------------------------------------------------------------------------------------
#
#                  p i n d e l ,   s t r a n d   c o u n t s
#
@test "pindel -sc" {
  run $PINDEL
  assert_failure 1
  assert_output --partial 'report_strand_counts'
}


# --------------------------------------------------------------------------------------
#
#                  p i n d e l ,   d e m o,   s t r a n d   c o u n t s
#
@test "pindel -sc: demo BAMs ➔ D, TD, INV, SI" {
  run $PINDEL -i simulated_config.txt -f $DEMOREF -o result/strand-counts -c ALL -sc
  assert_success
}

@test "pindel -sc: demo D" {
  run cmp result/strand-counts_D expected/strand-counts_D
  assert_success
}

@test "pindel -sc: demo TD" {
  run cmp result/strand-counts_TD expected/strand-counts_TD
  assert_success
}

@test "pindel -sc: demo INV" {
  run cmp result/strand-counts_INV expected/strand-counts_INV
  assert_success
}

@test "pindel -sc: demo SI" {
  run cmp result/strand-counts_SI expected/strand-counts_SI
  assert_success
}


# --------------------------------------------------------------------------------------
#
#                  p i n d e l ,   r e a l,   s t r a n d   c o u n t s
#
@test "pindel -sc: real BAMs ➔ D, TD, INV, SI" {
  run $PINDEL -i real_config.txt -f $REALREF -o result/strand-counts-real -c ALL -sc
  assert_success
}

@test "pindel -sc: real D" {
  run cmp result/strand-counts-real_D expected/strand-counts-real_D
  assert_success
}

@test "pindel -sc: real INV" {
  run cmp result/strand-counts-real_INV expected/strand-counts-real_INV
  assert_success
}


# --------------------------------------------------------------------------------------
#
#                          p i n d e l 2 v c f ,   d e m o
#
@test "pindel2vcf: demo D ➔ D.vcf" {
  run $PINDEL2VCF -e 5 -r $DEMOREF -mc 1 -R HUMAN_G1K_V2 -d 20100101 -p input/no-strand-counts_D --vcf result/no-strand-counts_D.vcf
  assert_success
}

@test "pindel2vcf: demo D.vcf" {
  run cmp result/no-strand-counts_D.vcf expected/no-strand-counts_D.vcf
  assert_success
}

@test "pindel2vcf: demo TD ➔ demo TD.vcf" {
  run $PINDEL2VCF -e 5 -r $DEMOREF -mc 1 -R HUMAN_G1K_V2 -d 20100101 -p input/no-strand-counts_TD --vcf result/no-strand-counts_TD.vcf
  assert_success
}

@test "pindel2vcf: demo TD.vcf" {
  run cmp result/no-strand-counts_TD.vcf expected/no-strand-counts_TD.vcf
  assert_success
}

@test "pindel2vcf: demo INV ➔ demo INV.vcf" {
  run $PINDEL2VCF -e 5 -r $DEMOREF -mc 1 -R HUMAN_G1K_V2 -d 20100101 -p input/no-strand-counts_INV --vcf result/no-strand-counts_INV.vcf
  assert_success
}

@test "pindel2vcf: demo INV.vcf" {
  run cmp result/no-strand-counts_INV.vcf expected/no-strand-counts_INV.vcf
  assert_success
}

@test "pindel2vcf: demo SI ➔ demo SI.vcf" {
  run $PINDEL2VCF -e 5 -r $DEMOREF -mc 1 -R HUMAN_G1K_V2 -d 20100101 -p input/no-strand-counts_SI --vcf result/no-strand-counts_SI.vcf
  assert_success
}

@test "pindel2vcf: demo SI.vcf" {
  run cmp result/no-strand-counts_SI.vcf expected/no-strand-counts_SI.vcf
  assert_success
}


# --------------------------------------------------------------------------------------
#
#                          p i n d e l 2 v c f ,   r e a l
#
#
@test "pindel2vcf: real D ➔ D.vcf" {
  run $PINDEL2VCF -e 2 -r $REALREF -mc 1 -R HUMAN_G1K_V2 -d 20100101 -p input/no-strand-counts-real_D --vcf result/no-strand-counts-real_D.vcf
  assert_success
}

@test "pindel2vcf: real D.vcf" {
  run cmp result/no-strand-counts-real_D.vcf expected/no-strand-counts-real_D.vcf
  assert_success
}

@test "pindel2vcf: real INV ➔ INV.vcf" {
  run $PINDEL2VCF -e 2 -r $REALREF -mc 1 -R HUMAN_G1K_V2 -d 20100101 -p input/no-strand-counts-real_INV --vcf result/no-strand-counts-real_INV.vcf
  assert_success
}

@test "pindel2vcf: real INV.vcf" {
  run cmp result/no-strand-counts-real_INV.vcf expected/no-strand-counts-real_INV.vcf
  assert_success
}


# --------------------------------------------------------------------------------------
#
#                   p i n d e l 2 v c f ,   s t r a n d   c o u n t s
#
#
@test "pindel2vcf -sc" {
  run $PINDEL2VCF
  assert_success
  assert_output --partial '-sc'
}

# --------------------------------------------------------------------------------------
#
#             p i n d e l 2 v c f ,   d e m o,   s t r a n d   c o u n t s
#
#
@test "pindel2vcf -sc: demo D ➔ D.vcf" {
  run $PINDEL2VCF -e 5 -r $DEMOREF -mc 1 -sc -R HUMAN_G1K_V2 -d 20100101 -p input/strand-counts_D --vcf result/strand-counts_D.vcf
  assert_success
}

@test "pindel2vcf -sc: demo D.vcf" {
  run cmp result/strand-counts_D.vcf expected/strand-counts_D.vcf
  assert_success
}

@test "pindel2vcf -sc: demo TD ➔ TD.vcf" {
  run $PINDEL2VCF -e 5 -r $DEMOREF -mc 1 -sc -R HUMAN_G1K_V2 -d 20100101 -p input/strand-counts_TD --vcf result/strand-counts_TD.vcf
  assert_success
}

@test "pindel2vcf -sc: demo TD.vcf" {
  run cmp result/strand-counts_TD.vcf expected/strand-counts_TD.vcf
  assert_success
}

@test "pindel2vcf -sc: demo INV ➔ INV.vcf" {
  run $PINDEL2VCF -e 5 -r $DEMOREF -mc 1 -sc -R HUMAN_G1K_V2 -d 20100101 -p input/strand-counts_INV --vcf result/strand-counts_INV.vcf
  assert_success
}

@test "pindel2vcf -sc: demo INV.vcf" {
  run cmp result/strand-counts_INV.vcf expected/strand-counts_INV.vcf
  assert_success
}

@test "pindel2vcf -sc: demo SI ➔ SI.vcf" {
  run $PINDEL2VCF -e 5 -r $DEMOREF -mc 1 -sc -R HUMAN_G1K_V2 -d 20100101 -p input/strand-counts_SI --vcf result/strand-counts_SI.vcf
  assert_success
}

@test "pindel2vcf -sc: demo SI.vcf" {
  run cmp result/strand-counts_SI.vcf expected/strand-counts_SI.vcf
  assert_success
}


# --------------------------------------------------------------------------------------
#
#             p i n d e l 2 v c f ,   r e a l,   s t r a n d   c o u n t s
#
#
@test "pindel2vcf -sc: real D ➔ D.vcf" {
  run $PINDEL2VCF -e 2 -r $REALREF -mc 1 -sc -R HUMAN_G1K_V2 -d 20100101 -p input/strand-counts-real_D --vcf result/strand-counts-real_D.vcf
  assert_success
}

@test "pindel2vcf -sc: real D.vcf" {
  run cmp result/strand-counts-real_D.vcf expected/strand-counts-real_D.vcf
  assert_success
}

@test "pindel2vcf -sc: real INV ➔ INV.vcf" {
  run $PINDEL2VCF -e 2 -r $REALREF -mc 1 -sc -R HUMAN_G1K_V2 -d 20100101 -p input/strand-counts-real_INV --vcf result/strand-counts-real_INV.vcf
  assert_success
}

@test "pindel2vcf -sc: real INV.vcf" {
  run cmp result/strand-counts-real_INV.vcf expected/strand-counts-real_INV.vcf
  assert_success
}


