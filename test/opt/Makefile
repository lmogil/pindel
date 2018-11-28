.PHONY: test

define TESTREF
hs_ref_chr5.fa
endef

define TESTREFINDEX
hs_ref_chr5.fa.fai
endef

define nl


endef

all: test remove-empties

ifeq ($(wildcard $(TESTREF)),)
  $(error Please provide an indexed reference FASTA named '$(TESTREF)',$(nl)                 containing complete human Chromosome 5 named '5')
endif

ifeq ($(wildcard $(TESTREFINDEX)),)
  $(error Please index the test reference FASTA by running 'samtools faidx $(TESTREF)')
endif

test: hs_ref_chr5.fa hs_ref_chr5.fa.fai
	mkdir -p result
	@./suite.bats

remove-empties:
	@rm -f \
    result/no-strand-counts-real_BP \
    result/no-strand-counts-real_CloseEndMapped \
    result/no-strand-counts-real_INT_final \
    result/no-strand-counts-real_LI \
    result/no-strand-counts-real_RP \
    result/no-strand-counts-real_SI \
    result/no-strand-counts-real_TD \
    result/no-strand-counts_BP \
    result/no-strand-counts_CloseEndMapped \
    result/no-strand-counts_INT_final \
    result/no-strand-counts_LI \
    result/no-strand-counts_RP \
    result/strand-counts-real_BP \
    result/strand-counts-real_CloseEndMapped \
    result/strand-counts-real_INT_final \
    result/strand-counts-real_LI \
    result/strand-counts-real_RP \
    result/strand-counts-real_SI \
    result/strand-counts-real_TD \
    result/strand-counts_BP \
    result/strand-counts_CloseEndMapped \
    result/strand-counts_INT_final \
    result/strand-counts_LI \
    result/strand-counts_RP

clean:
	@rm -f result/*
