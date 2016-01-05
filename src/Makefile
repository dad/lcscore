PY = $$PYTHONEXE
GENOMES = ../../../../Dropbox/bio/genomes
DEV = ../../../../develop
DATA = ../data
SGAGG_DATA = ../../sg-agg/data
FIG = ../figures
TMP = ../../../../tmp
LIB = ../../lib/src


#PDOMAIN_START = 880
#PDOMAIN_END = 1450
#RRMS_START = 80
#RRMS_END = 879
PDOMAIN_START = 780
PDOMAIN_END = 1206
RRMS_START = 80
RRMS_END = 779

NSP_FG_START = 1
NSP_FG_END = 1510
NSP_STRUCT_START = 1511
NSP_STRUCT_END = 2126

PBP1_LSMAD_START = 640 
PBP1_LSMAD_END = 825
PBP1_PB_START = 1591
PBP1_PB_END = 2160

all: run

hmm:
	$(PY) hmm_test.py

strevol-test:
	$(PY) structevol_test.py

strevol:
	$(PY) structure-from-alignment.py $(SGAGG_DATA)/alignment/pab1-all.fa --out $(DATA)/pab1-all-profile.txt
	$(PY) structure-from-alignment.py $(SGAGG_DATA)/alignment/pbp1-all-aligned.fa --out $(DATA)/pbp1-all-aligned-profile.txt
	$(PY) structure-from-alignment.py $(SGAGG_DATA)/alignment/nsp1-all-aligned.fa --out $(DATA)/nsp1-all-aligned-profile.txt
	$(PY) structure-from-alignment.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --out $(DATA)/nup57-fungal-final-profile.txt

# fungal eIF4G
# fungal importin-beta
ORTHOID = EOG7F7WJ5
ORTHONAME = kap123-fungal
ORTHOID = EOG7SFJ6V
ORTHONAME = fip1-fungal
ORTHOID = EOG7NKKW1
ORTHONAME = cyc8-fungal
ORTHOID = EOG7MKWFN
ORTHONAME = gus1-fungal
ORTHOID = EOG7CG77Z
ORTHONAME = mes1-fungal
ORTHOID = EOG7RBZK2
ORTHONAME = tef4-fungal
ORTHOID = EOG71ZPBR
ORTHONAME = arc1-fungal
ORTHOID =  EOG7TJ3T2
ORTHONAME = ded1-fungal
ORTHOANCHOR = DED1_YEAST
ORTHOID = EOG7ZPNWZ
ORTHONAME = nup57-fungal
ORTHOANCHOR = YGR119C
ORTHOID = EOG728916
ORTHONAME = ssa1-fungal
ORTHOANCHOR = SSA1
ORTHOID = EOG7PCJGF
ORTHONAME = hsc70-metazoa
ORTHOANCHOR = P11142
ORTHOID = EOG6GR442
ORTHONAME = dnak-bacteria
ORTHOANCHOR = DNAK_ECOLI

ORTHOID = EOG7PZRZB
ORTHONAME = p62-nup-metazoan
ORTHOANCHOR = P37198
# hnRNP-G
ORTHOID = EOG780RPD
ORTHONAME = hnrnpg-metazoan
ORTHOANCHOR = P38159


TMPID = tmpo
ORTHOID = EOG706125
ORTHONAME = pub1-fungal
ORTHOANCHOR = YNL016W

ORTHOID = EOG7SFJ7D
ORTHONAME = sis1-fungal
ORTHOANCHOR = SIS1

ORTHOID = EOG77M90H
ORTHONAME = ydj1-fungal
ORTHOANCHOR = P25491

 
# metazoan eIF4G
#ORTHOID = EOG7D59N2
#ORTHONAME = eif4g
# fungal eIF4G
#ORTHOID = EOG7VTDWG
#ORTHONAME = eif4g-fungal

orthodb: odbfetch odbfilter odbalign odbprune

odbfetch:
	$(PY) fetch-orthodb.py $(ORTHOID) --fasta-out $(SGAGG_DATA)/alignment/$(TMPID)-orthodb.fa
odbfilter:
	$(PY) filter-orthodb.py $(SGAGG_DATA)/alignment/$(TMPID)-orthodb.fa --orthodb-names ../../data/OrthoDB7_SPECIES_tabtext --sentinel $(ORTHOANCHOR) --remove-x --out $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb-filtered.log --fasta-out $(SGAGG_DATA)/alignment/$(TMPID)-orthodb-filtered.fa
odbalign:
	$(PY) $(LIB)/muscle.py $(SGAGG_DATA)/alignment/$(TMPID)-orthodb-filtered.fa --out $(SGAGG_DATA)/alignment/$(TMPID)-orthodb-filtered-aligned.fa
odbprune:
	$(PY) prune-alignment.py $(SGAGG_DATA)/alignment/$(TMPID)-orthodb-filtered-aligned.fa --anchor $(ORTHOANCHOR) -1 --coverage 0.3 --identity 0.2 --degap --out $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb-aligned-pruned.log --fasta-out $(SGAGG_DATA)/alignment/$(ORTHONAME)-final.fa

# 
# $(PY) prune-alignment.py $(SGAGG_DATA)/alignment/pab1-total.fa --anchor YER165W -1 --degap --fasta-out $(SGAGG_DATA)/alignment/pab1-all.fa

fetch:
#	$(PY) fetch-orthodb.py --orthodb-id EOG7D59X8 --name nab3 --out-dir $(SGAGG_DATA)/alignment/
	$(PY) fetch-orthodb.py $(ORTHOID) --fasta-out $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb.fa
#	$(PY) fetch-orthodb.py EOG7RV9FP --fasta-out $(SGAGG_DATA)/alignment/pab1vtest-orthodb.fa

filter:
#	$(PY) filter-orthodb.py $(SGAGG_DATA)/alignment/hsp70-orthodb-metazoa.fa ../../data/OrthoDB7_SPECIES_tabtext $(DATA)/preferred-strains.txt --include cognate --sequence EEVD --out $(SGAGG_DATA)/alignment/hsp70-metazoa.log --fasta-out $(SGAGG_DATA)/alignment/hsp70-metazoa.fa
#	$(PY) filter-orthodb.py $(SGAGG_DATA)/alignment/npl3-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext $(DATA)/preferred-strains.txt --out $(SGAGG_DATA)/alignment/npl3-orthodb.log --fasta-out $(SGAGG_DATA)/alignment/npl3-orthodb-filtered.fa
#	$(PY) filter-orthodb.py $(DATA)/allpab.fa ../../data/OrthoDB7_SPECIES_tabtext $(DATA)/preferred-strains.txt --include IPR002004 --include PABP --include "polyadenylate-binding protein 1" --min-length 400 --exclude splicing --out $(SGAGG_DATA)/alignment/pab1test-orthodb.log --fasta-out $(SGAGG_DATA)/alignment/pab1test-orthodb-filtered.fa
#	$(PY) filter-orthodb.py $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext $(DATA)/preferred-strains.txt --remove-x --out $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb.log --fasta-out $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb-filtered.fa
	$(PY) filter-orthodb.py $(DATA)/allpab.fa ../../data/OrthoDB7_SPECIES_tabtext --sentinel PABP_YEAST --include IPR002004 --include PABP --include "polyadenylate-binding protein 1" --min-length 400 --out $(SGAGG_DATA)/alignment/pab1test-orthodb.log --fasta-out $(DATA)/pab1test-orthodb-filtered.fa

align:
#	$(PY) $(LIB)/muscle.py $(DATA)/pab1test-orthodb-filtered.fa --out $(DATA)/pab1test-orthodb-filtered-aligned.fa
	$(PY) $(LIB)/muscle.py $(DATA)/pab1test-orthodb-filtered.fa --out $(DATA)/pab1test-orthodb-filtered-aligned.fa

rename:
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/p62-nup-metazoan-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext --align --remove-x --out $(SGAGG_DATA)/alignment/p62-nup-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/nsp1-fungal-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext --align --remove-x --out $(SGAGG_DATA)/alignment/nsp1-fungal-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/nsp1-fungal-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext -f $(SGAGG_DATA)/alignment/p62-nup-metazoan-orthodb.fa --align --remove-x --out $(SGAGG_DATA)/alignment/nsp1-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/nsp1-fungal-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext -f $(SGAGG_DATA)/alignment/p62-nup-metazoan-orthodb.fa --align --remove-x --out $(SGAGG_DATA)/alignment/test.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/nsp1-test.fa ../../data/OrthoDB7_SPECIES_tabtext --align --remove-x --out $(SGAGG_DATA)/alignment/test.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/elastin-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext --align --remove-x --out $(SGAGG_DATA)/alignment/elastin-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/pbp1-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext --align --remove-x --out $(SGAGG_DATA)/alignment/pbp1-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/sro9-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext --align --remove-x --out $(SGAGG_DATA)/alignment/sro9-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/src-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext --align --remove-x --out $(SGAGG_DATA)/alignment/src-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/nab3-orthodb.fa ../../data/OrthoDB7_SPECIES_tabtext --align --remove-x --out $(SGAGG_DATA)/alignment/nab3-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/npl3-orthodb-filtered.fa ../../data/OrthoDB7_SPECIES_tabtext --align --out $(SGAGG_DATA)/alignment/npl3-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/pab1test-orthodb-filtered.fa ../../data/OrthoDB7_SPECIES_tabtext --align --out $(SGAGG_DATA)/alignment/pab1test-orthodb-aligned.fa
#	$(PY) rename-orthodb-hits.py $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb-filtered.fa ../../data/OrthoDB7_SPECIES_tabtext --align --out $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb-aligned.fa

prune:
#	$(PY) prune-species.py $(SGAGG_DATA)/alignment/nsp1-orthodb-aligned.fa $(DATA)/preferred-strains.txt --gap-thresh 0.2 --id-thresh 0.2 --coverage-thresh 0.5 --out $(SGAGG_DATA)/alignment/nsp1-all-aligned.log --fasta-out $(SGAGG_DATA)/alignment/nsp1-all-aligned.fa
#	$(PY) prune-species.py $(SGAGG_DATA)/alignment/elastin-orthodb-aligned.fa $(DATA)/preferred-strains.txt --gap-thresh 0.2 --id-thresh 0.2 --coverage-thresh 0.5 --out $(SGAGG_DATA)/alignment/elastin-all-aligned.log --fasta-out $(SGAGG_DATA)/alignment/elastin-all-aligned.fa
#	$(PY) prune-species.py $(SGAGG_DATA)/alignment/pbp1-orthodb-aligned.fa $(DATA)/preferred-strains.txt --gap-thresh 0.2 --id-thresh 0.2 --coverage-thresh 0.5 --out $(SGAGG_DATA)/alignment/pbp1-all-aligned.log --fasta-out $(SGAGG_DATA)/alignment/pbp1-all-aligned.fa
#	$(PY) prune-species.py $(SGAGG_DATA)/alignment/nab3-orthodb-aligned.fa $(DATA)/preferred-strains.txt --gap-thresh 0.2 --id-thresh 0.2 --coverage-thresh 0.5 --out $(SGAGG_DATA)/alignment/nab3-all-aligned.log --fasta-out $(SGAGG_DATA)/alignment/nab3-all-aligned.fa
#	$(PY) prune-species.py $(SGAGG_DATA)/alignment/hsp70-orthodb.fa $(DATA)/preferred-strains.txt --out $(SGAGG_DATA)/alignment/hsp70-all.log --fasta-out $(SGAGG_DATA)/alignment/hsp70-all.fa
#	$(PY) prune-species.py $(SGAGG_DATA)/alignment/npl3-orthodb-aligned.fa $(DATA)/preferred-strains.txt --out $(SGAGG_DATA)/alignment/npl3-aligned-pruned.log --fasta-out $(SGAGG_DATA)/alignment/npl3-aligned-pruned.fa
#	$(PY) prune-species.py $(SGAGG_DATA)/alignment/pab1test-orthodb-aligned.fa $(DATA)/preferred-strains.txt --out $(SGAGG_DATA)/alignment/pab1test-aligned-pruned.log --fasta-out $(SGAGG_DATA)/alignment/pab1test-aligned-pruned.fa
#	$(PY) prune-species.py $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb-aligned.fa $(DATA)/preferred-strains.txt --out $(SGAGG_DATA)/alignment/$(ORTHONAME)-aligned-pruned.log --fasta-out $(SGAGG_DATA)/alignment/$(ORTHONAME)-aligned-pruned.fa
#	$(PY) prune-alignment.py $(DATA)/pab1test-orthodb-filtered-aligned.fa --anchor PABP_YEAST --coverage 0.9 --identity 0.4 --degap --out $(DATA)/pab1-orthodb-aligned-pruned.log --fasta-out $(DATA)/pab1-orthodb-aligned-pruned.fa
	$(PY) prune-alignment.py $(SGAGG_DATA)/alignment/$(ORTHONAME)-orthodb-aligned.fa --anchor PABP_YEAST --coverage 0.9 --identity 0.4 --degap --out $(DATA)/pab1-orthodb-aligned-pruned.log --fasta-out $(DATA)/pab1-orthodb-aligned-pruned.fa


pmutants:
#	$(PY) proteins-to-genes.py $(DATA)/pab1p-mutants.fa --suggest $(DATA)/scer-pab1p-suggestions.fa --prefix TCAATTGGCTCAACAAATCCAAGCCAGAAATCAAATGAGA --suffix CAAAAGCAAAGACAAGCTTTGGGTGAACAATTATACAAGA --out $(DATA)/pab1p-mutant-gblocks.txt
#	$(PY) proteins-to-genes.py $(DATA)/pab1p-mutants.fa --suggest $(DATA)/scer-pab1p-suggestions-avoid.fa --prefix TCAATTGGCTCAACAAATCCAAGCCAGAAATCAAATGAGA --suffix CAAAAGCAAAGACAAGCTTTGGGTGAACAATTATACAAGA --out $(DATA)/pab1p-mutant-gblocks-avoid.txt
#	$(PY) proteins-to-genes.py $(DATA)/pab1p-hydrophobic-mutants.fa --suggest $(DATA)/scer-pab1p-suggestions.fa --prefix TCAATTGGCTCAACAAATCCAAGCCAGAAATCAAATGAGA --suffix CAAAAGCAAAGACAAGCTTTGGGTGAACAATTATACAAGA --out $(DATA)/pab1p-hydro-mutant-gblocks.txt
	$(PY) proteins-to-genes.py $(DATA)/pab1p-mutants-v2.fa --suggest $(DATA)/scer-pab1p-suggestions-avoid.fa --prefix ATGAGATATCAACAAGCCACGGCCGCAGCAGCAGCAGCGGCAGCA --suffix AGGAATGCAAATGACAATAATCAATTTTATCAACAAAAGCAAAGGC --out $(DATA)/pab1p-mutants-v2.txt

droplet:
	$(PY) droplet_driver.py $(SGAGG_DATA)/alignment/pab1-all.fa --packing-density 0.72 --degap --sentinel-species S.cerevisiae --start-sequence QLAQQIQ --end-sequence ANDNNQFYQ --aas FY -a KR --aas ILMV -a AGNPQ -a DE -a ST -a HCW --out $(DATA)/pab1p-droplet-structure.txt
	$(PY) droplet_driver.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --packing-density 0.72 --degap --sentinel-species S.cerevisiae --start-sequence FSFGTANNNS --end-sequence DSLKLNSKPVE --aas FY -a KR --aas ILMV -a AGNPQ -a DE -a ST -a HCW --out $(DATA)/nsp1-fg-droplet-structure.txt
	$(PY) droplet_driver.py $(SGAGG_DATA)/alignment/pbp1-all-aligned.fa --packing-density 0.72 --degap --sentinel-species S.cerevisiae --start-sequence QQRQLNSMGNAV  --end-sequence MMGSGPSGMP --aas FY -a KR --aas ILMV -a AGNPQ -a DE -a ST -a HCW --out $(DATA)/pbp1-pb-droplet-structure.txt
	$(PY) droplet_driver.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --packing-density 0.72 --degap --sentinel-species S.cerevisiae --start-sequence MFGFSGSNNGFGNKPAG --end-sequence QQQQLSNYPQQIQEQ --aas FY -a KR --aas ILMV -a AGNPQ -a DE -a ST -a HCW --out $(DATA)/nup57-fg-droplet-structure.txt

sim:
#	$(PY) aggsim-driver.py 10 45 --monomer-radius 1.0 \
#	--k-nucleation 1e-13 --k-denucleation 1e-6 \
#	--k-assembly 1e-8 -T 42 --output-step 100 --out $(DATA)/sim-test.txt
	$(PY) aggsim-driver.py 10 45 --monomer-radius 1.0 \
	--k-nucleation 1e-15 --k-denucleation 1e-3 \
	--k-assembly 1e-8 --k-disassembly 1e-6 -T 30 --jump-time 1000 --jump-temp 42 --output-time 60 --output-step 100 --out $(DATA)/sim-nucleation-1e-15.txt
#	$(PY) aggsim-driver.py 10 45 --monomer-radius 1.0 \
#	--k-nucleation 1e-14 --k-denucleation 1e-3 \
#	--k-assembly 1e-8 --k-disassembly 1e-6 -T 30 --jump-time 1000 --jump-temp 42 --output-time 60 --output-step 100 --out $(DATA)/sim-nucleation-1e-14.txt
#	$(PY) aggsim-driver.py 10 45 --monomer-radius 1.0 \
#	--k-nucleation 1e-12 --k-denucleation 1e-3 \
#	--k-assembly 1e-8 --k-disassembly 1e-6 -T 30 --jump-time 1000 --jump-temp 42 --output-time 60 --output-step 100 --out $(DATA)/sim-nucleation-1e-12.txt
#	$(PY) aggsim-driver.py 10 45 --monomer-radius 4.5 \
#	--k-nucleation 1e-13 --k-denucleation 1e-3 \
#	--k-assembly 1e-8 --k-disassembly 1e-6 -T 30 --jump-time 2000 --jump-temp 40 --output-step 10 --out $(DATA)/sim-test.txt

simtest:
	$(PY) aggsim_test.py

MOTLEN = 6
findmax:
#	$(PY) find-max-composition.py $(GENOMES)/scer_orf_trans.fasta 10 --start 406 --end 502 --query-gene PAB1
#	$(PY) find-max-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 10 --degap --query-orf S.cerevisiae --start 567 --end 879
#	$(PY) find-max-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 10 --degap --query-orf H.sapiens --start 567 --end 879
#	$(PY) find-max-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 10 --degap --query-orf D.melanogaster --start 567 --end 879
#	$(PY) find-max-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 10 --degap --query-orf S.pombe --start 567 --end 879
	$(PY) find-max-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 6  --degap --start 650 --end 1090 --out $(DATA)/pab1p-max-composition-L06.txt
	$(PY) find-max-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 8  --degap --start 650 --end 1090 --out $(DATA)/pab1p-max-composition-L08.txt
	$(PY) find-max-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 10 --degap --start 650 --end 1090 --out $(DATA)/pab1p-max-composition-L10.txt
	$(PY) find-max-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 12 --degap --start 650 --end 1090 --out $(DATA)/pab1p-max-composition-L12.txt
	$(PY) find-max-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 14 --degap --start 650 --end 1090 --out $(DATA)/pab1p-max-composition-L14.txt
# --out $(DATA)/test.txt

findmatch:
	$(PY) find-by-composition.py $(GENOMES)/scer_orf_trans.fasta 50 --threshold 0.85 --motif QLAQQIQARNQMRYQQATAAAAAAAAGMPGQFMPPMFYGVMPPRGVPFNGPNPQQMNPMGGMPKNGMPPQFRNGPVYGVPPQGGFPRNANDNNQFYQ --out $(DATA)/scer-pab1p-matches-50.txt
#	$(PY) find-by-composition.py $(GENOMES)/scer_orf_trans.fasta 20 --threshold 0.85 --motif QLAQQIQARNQMRYQQATAAAAAAAAGMPGQFMPPMFYGVMPPRGVPFNGPNPQQMNPMGGMPKNGMPPQFRNGPVYGVPPQGGFPRNANDNNQFYQ --out $(DATA)/scer-pab1p-matches.txt
#	$(PY) find-by-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 20 --threshold 0.7 --motif AAAGMPPQRT --degap --out $(DATA)/all-species-AAAGMPPQRT-matches-test.txt
#	$(PY) find-by-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 20 --threshold 0.7 --motif FGGMNPPPQR --degap --out $(DATA)/all-species-FGGMNPPPQR-matches-test.txt
#	$(PY) find-by-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa 20 --threshold 0.7 --motif NAGRGMP --degap --out $(DATA)/all-species-NAGRGMP-matches-test.txt

comp:
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence QLAQQIQ --end-sequence ANDNNQFYQ -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/pab1p-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence SASLYVGDLEPS --end-sequence AIAQRKDVRRS -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/pab1rrm-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pbp1-all-aligned.fa --degap --sentinel-species S.cerevisiae --start-sequence QQRQLNSMGNAV  --end-sequence MMGSGPSGMP -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/pbp1-pb-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pbp1-all-aligned.fa --degap --sentinel-species S.cerevisiae --start-sequence DDRQDYLLANSI --end-sequence DEEDLYSGVDRR -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/pbp1-lsmad-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence FSFGTANNNS --end-sequence DSLKLNSKPVE -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/nsp1-fg-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence TLDDLVTKWTNQL --end-sequence SLEKQINSI -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/nsp1-struct-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence MFGFSGSNNGFGNKPAG --end-sequence QQQQLSNYPQQIQEQ -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/nup57-fg-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence YLNEVLEKDAAIVK --end-sequence YLNEVLEKDAAIVK -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/nup57-struct-composition.txt

sis-comp:
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/sis1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence GGPSFGPGGPGG --end-sequence FSNEDAFNIFSQ -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas FG --aas GM --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/sis1-fg-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/sis1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence FFGGSSPFGGADDS --end-sequence SPTYPEEE -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas FG --aas GM --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/sis1-gm-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/sis1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence GGPSFGPGGPGG --end-sequence SPTYPEEE -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas FG --aas GM --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/sis1-lc-composition.txt
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/sis1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence TVQVNLPVSLEDLF --end-sequence PISLNDAQKRAIDENF -s ACDEFGHIKLMNPQRSTVWY \
		--aas FY -a KR --aas FG --aas GM --aas AGNPQ -a ILMVW -a ST -a DE --out $(DATA)/sis1-struct-composition.txt

motif:
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --aas FY --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY-composition-motif.txt --isolate-out $(DATA)/pab1p-FY-motifs.fa
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --aas FY --aas KR --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_KR-composition-motif.txt --isolate-out $(DATA)/pab1p-FY_KR-motifs.fa
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --aas FY --aas P --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_P-composition-motif.txt --isolate-out $(DATA)/pab1p-FY_P-motifs.fa
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --aas FY --aas QN --start $(PDOMAIN_START) --end $(PDOMAIN_END) --isolate-out $(DATA)/pab1p-FY_QN-motifs.fa
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --aas FY --aas PQ --start $(PDOMAIN_START) --end $(PDOMAIN_END) --isolate-out $(DATA)/pab1p-FY_PQ-motifs.fa
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --aas FY --aas Q --start $(PDOMAIN_START) --end $(PDOMAIN_END) --isolate-out $(DATA)/pab1p-FY_Q-motifs.fa
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --aas FY --aas ILMVW --start $(PDOMAIN_START) --end $(PDOMAIN_END) --isolate-out $(DATA)/pab1p-FY_ILMVW-motifs.fa
	$(PY) quantify-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --aas FY --aas G --start $(PDOMAIN_START) --end $(PDOMAIN_END) --isolate-out $(DATA)/pab1p-FY_G-motifs.fa


sep:
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY Q --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_Q-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY P --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_P-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY PQ --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_PQ-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY PQK --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_PQK-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY PQKS --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_PQKS-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY PQKNESD --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_PQKNESD-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY PQKNESDG --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_PQKNESDG-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY ACDEGHIKLMNPQRSTVW --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_ACDEGHIKLMNPQRSTVW-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY M --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_M-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY G --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_G-separator-distributions.txt
	$(PY) quantify-separators.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa FY N --max-dist 20 --degap --start $(PDOMAIN_START) --end $(PDOMAIN_END) --out $(DATA)/pab1p-FY_N-separator-distributions.txt

dist: dist-pab dist-nsp dist-nup

dist-pab:
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence QLAQQIQ --end-sequence ANDNNQFYQ --aas FY --out $(DATA)/pab1p-FY-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence QLAQQIQ --end-sequence ANDNNQFYQ --aas GKR --out $(DATA)/pab1p-GKR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence QLAQQIQ --end-sequence ANDNNQFYQ --aas ST --out $(DATA)/pab1p-ST-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence QLAQQIQ --end-sequence ANDNNQFYQ --aas KR --out $(DATA)/pab1p-KR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence QLAQQIQ --end-sequence ANDNNQFYQ --aas ILMVW --out $(DATA)/pab1p-ILMVW-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence QLAQQIQ --end-sequence ANDNNQFYQ --aas AGNPQ --out $(DATA)/pab1p-AGNPQ-distances.txt

	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence SASLYVGDLEPS --end-sequence AIAQRKDVRRS --aas FY --out $(DATA)/pab1rrm-FY-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence SASLYVGDLEPS --end-sequence AIAQRKDVRRS --aas GKR --out $(DATA)/pab1rrm-GKR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence SASLYVGDLEPS --end-sequence AIAQRKDVRRS --aas ST --out $(DATA)/pab1rrm-ST-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence SASLYVGDLEPS --end-sequence AIAQRKDVRRS --aas KR --out $(DATA)/pab1rrm-KR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence SASLYVGDLEPS --end-sequence AIAQRKDVRRS --aas ILMVW --out $(DATA)/pab1rrm-ILMVW-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/pab1-all.fa --degap --sentinel-species S.cerevisiae --start-sequence SASLYVGDLEPS --end-sequence AIAQRKDVRRS --aas AGNPQ --out $(DATA)/pab1rrm-AGNPQ-distances.txt

dist-nsp:
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence FSFGTANNNS --end-sequence DSLKLNSKPVE --aas FY --out $(DATA)/nsp1-fgrepeat-FY-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence FSFGTANNNS --end-sequence DSLKLNSKPVE --aas FYKR --out $(DATA)/nsp1-fgrepeat-FYKR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence FSFGTANNNS --end-sequence DSLKLNSKPVE --aas ST --out $(DATA)/nsp1-fgrepeat-ST-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence FSFGTANNNS --end-sequence DSLKLNSKPVE --aas DE --out $(DATA)/nsp1-fgrepeat-DE-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence FSFGTANNNS --end-sequence DSLKLNSKPVE --aas KR --out $(DATA)/nsp1-fgrepeat-KR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence FSFGTANNNS --end-sequence DSLKLNSKPVE --aas ILMVW --out $(DATA)/nsp1-fgrepeat-ILMVW-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence FSFGTANNNS --end-sequence DSLKLNSKPVE --aas AGNPQ --out $(DATA)/nsp1-fgrepeat-AGNPQ-distances.txt

	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence TLDDLVTKWTNQL --end-sequence SLEKQINSI --aas FY --out $(DATA)/nsp1-struct-FY-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence TLDDLVTKWTNQL --end-sequence SLEKQINSI --aas FYKR --out $(DATA)/nsp1-struct-FYKR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence TLDDLVTKWTNQL --end-sequence SLEKQINSI --aas ST --out $(DATA)/nsp1-struct-ST-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence TLDDLVTKWTNQL --end-sequence SLEKQINSI --aas DE --out $(DATA)/nsp1-struct-DE-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence TLDDLVTKWTNQL --end-sequence SLEKQINSI --aas KR --out $(DATA)/nsp1-struct-KR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence TLDDLVTKWTNQL --end-sequence SLEKQINSI --aas ILMVW --out $(DATA)/nsp1-struct-ILMVW-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nsp1-p62-all.fa --degap --sentinel-species S.cerevisiae --start-sequence TLDDLVTKWTNQL --end-sequence SLEKQINSI --aas AGNPQ --out $(DATA)/nsp1-struct-AGNPQ-distances.txt

dist-nup:
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence MFGFSGSNNGFGNKPAG --end-sequence QQQQLSNYPQQIQEQ --aas FY --out $(DATA)/nup57-fgrepeat-FY-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence MFGFSGSNNGFGNKPAG --end-sequence QQQQLSNYPQQIQEQ --aas FYKR --out $(DATA)/nup57-fgrepeat-FYKR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence MFGFSGSNNGFGNKPAG --end-sequence QQQQLSNYPQQIQEQ --aas ST --out $(DATA)/nup57-fgrepeat-ST-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence MFGFSGSNNGFGNKPAG --end-sequence QQQQLSNYPQQIQEQ --aas DE --out $(DATA)/nup57-fgrepeat-DE-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence MFGFSGSNNGFGNKPAG --end-sequence QQQQLSNYPQQIQEQ --aas KR --out $(DATA)/nup57-fgrepeat-KR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence MFGFSGSNNGFGNKPAG --end-sequence QQQQLSNYPQQIQEQ --aas ILMVW --out $(DATA)/nup57-fgrepeat-ILMVW-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence MFGFSGSNNGFGNKPAG --end-sequence QQQQLSNYPQQIQEQ --aas AGNPQ --out $(DATA)/nup57-fgrepeat-AGNPQ-distances.txt

	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence VLKCKESWDPNTTKTKL --end-sequence NEVLEKDAAIVKK --aas FY --out $(DATA)/nup57-struct-FY-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence VLKCKESWDPNTTKTKL --end-sequence NEVLEKDAAIVKK --aas FYKR --out $(DATA)/nup57-struct-FYKR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence VLKCKESWDPNTTKTKL --end-sequence NEVLEKDAAIVKK --aas ST --out $(DATA)/nup57-struct-ST-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence VLKCKESWDPNTTKTKL --end-sequence NEVLEKDAAIVKK --aas DE --out $(DATA)/nup57-struct-DE-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence VLKCKESWDPNTTKTKL --end-sequence NEVLEKDAAIVKK --aas KR --out $(DATA)/nup57-struct-KR-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence VLKCKESWDPNTTKTKL --end-sequence NEVLEKDAAIVKK --aas ILMVW --out $(DATA)/nup57-struct-ILMVW-distances.txt
	$(PY) quantify-distances.py $(SGAGG_DATA)/alignment/nup57-fungal-final.fa --degap --sentinel-species S.cerevisiae --start-sequence VLKCKESWDPNTTKTKL --end-sequence NEVLEKDAAIVKK --aas AGNPQ --out $(DATA)/nup57-struct-AGNPQ-distances.txt


hsp70:
	$(PY) ../../../hsp70-binding/analysis/motif-finder.py --database $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --pssm rudiger --threshold 5.0 --report --out $(DATA)/pab1-motifs-r5.txt
	$(PY) ../../../hsp70-binding/analysis/motif-finder.py --database $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --pssm vandurme --threshold 5.0 --mask --out $(DATA)/pab1-motifs-v5.fa
	$(PY) ../../../hsp70-binding/analysis/motif-finder.py --database $(SGAGG_DATA)/alignment/pab1-orthologs.fa --degap --pssm vandurme --threshold 3.0 --mask --out $(DATA)/pab1-motifs-v3.fa
slidetest:
	$(PY) slide_test.py

MOT = SGNIFIKNLHPDIDNKALYDTFSVFGDILSSKIATDENGKSKGFGFVHFEEEGAAKEAIDALNGMLLNGQEIYVAPHLS
MOT = SSASLYVGDLEPSVSEAHLYDIFSPIGSVSSIRVCRDAITKTSLGYAYVNFNDHEAGRKAIEQLNYTPIKGRLCRIMWSQR
MOT = QLAQQIQARNQMRYQQATAAAAAAAAGMPGQFMPPMFYGVMPPRGVPFNGPNPQQMNPMGGMPKNGMPPQFRNGPVYGVPPQGGFPRNANDNNQFYQ
MOT = GMPGQFMPPMFYGVMPPRGVPFNGPNPQQMNPMGGMPKNGMPPQFRNGPVYGVPPQGGFPRNANDNNQFYQ
#
slide:
	$(PY) slide-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --motif $(MOT) --degap --query-orf S.cerevisiae --out $(DATA)/pab1-scer-slidecomp.txt
	$(PY) slide-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --motif $(MOT) --degap --query-orf O.parapolymorpha --out $(DATA)/pab1-opar-slidecomp.txt
	$(PY) slide-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --motif $(MOT) --degap --query-orf H.sapiens --out $(DATA)/pab1-hsap-slidecomp.txt
	$(PY) slide-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --motif $(MOT) --degap --query-orf X.laevis --out $(DATA)/pab1-xlae-slidecomp.txt
	$(PY) slide-composition.py $(SGAGG_DATA)/alignment/pab1-orthologs.fa --motif $(MOT) --degap --query-orf D.melanogaster --out $(DATA)/pab1-dmel-slidecomp.txt

run:
#	$(PY) gelscore_driver.py $(GENOMES)/scer_orf_trans.fasta --min-size 10 --threshold 0.9 --out $(DATA)/scer-regions.txt
	$(PY) gelscore_driver.py $(GENOMES)/dmel.cds.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/dmel-regions.txt &
	$(PY) gelscore_driver.py $(GENOMES)/dsec.cds.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/dsec-regions.txt &
	$(PY) gelscore_driver.py $(GENOMES)/dere.cds.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/dere-regions.txt &
	$(PY) gelscore_driver.py $(GENOMES)/dyak.cds.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/dyak-regions.txt &
	$(PY) gelscore_driver.py $(GENOMES)/dvir.cds.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/dvir-regions.txt &
#	$(PY) gelscore_driver.py $(GENOMES)/Sbay.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/sbay-regions.txt &
#	$(PY) gelscore_driver.py $(GENOMES)/Calb.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/calb-regions.txt &
#	$(PY) gelscore_driver.py $(GENOMES)/Cgla.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/cgla-regions.txt &
#	$(PY) gelscore_driver.py $(GENOMES)/Kwal.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/kwal-regions.txt &
#	$(PY) gelscore_driver.py $(GENOMES)/Scer.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/scer-regions2.txt &
#	$(PY) gelscore_driver.py $(GENOMES)/Spar.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/spar-regions.txt &
#	$(PY) gelscore_driver.py $(GENOMES)/Spom.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/spom-regions.txt &
#	$(PY) gelscore_driver.py $(GENOMES)/Smik.fasta --translate --min-size 10 --threshold 0.9 --out $(DATA)/smik-regions.txt &

runtest:
	$(PY) gelscore_driver.py $(GENOMES)/scer_orf_trans.fasta --min-size 10 --threshold 0.9 --debug --out $(DATA)/scer-regions-test.txt

test:
	$(PY) gelscore_test.py

TESTTHRESH = 0.025
TESTSIZE = 15
drivertest:
#	$(PY) gelscore_driver.py $(GENOMES)/scer_orf_trans.fasta -g --min-size $(TESTSIZE) --threshold 0.1 --out $(DATA)/scer-regions-test.txt
#	$(PY) gelscore_driver.py $(GENOMES)/scer_orf_trans.fasta --query-gene PAB1 --min-size $(TESTSIZE) --threshold $(TESTTHRESH)
#	$(PY) gelscore_driver.py $(GENOMES)/scer_orf_trans.fasta --query-gene PAB1 --query-gene PBP1 --query-gene SUP35 --query-gene RNQ1 \
#		--query-gene NEW1 --query-gene URE2 --query-gene NPL3 --query-gene YRA1 --query-gene NAB2 \
#		--min-size $(TESTSIZE) --threshold $(TESTTHRESH) --out $(DATA)/scer-region-subset.txt
	$(PY) gelscore_driver.py $(DATA)/ddx4.fa --min-size $(TESTSIZE) --threshold $(TESTTHRESH) --out $(DATA)/ddx4-scores.txt
#	$(PY) gelscore_driver.py $(DATA)/fus.fa --min-size $(TESTSIZE) --threshold $(TESTTHRESH) 
#	$(PY) gelscore_driver.py $(GENOMES)/scer_orf_trans.fasta --query-gene NPL3 --min-size $(TESTSIZE) --threshold $(TESTTHRESH) 
#	$(PY) gelscore_driver.py $(DATA)/taf15.fa --min-size $(TESTSIZE) --threshold $(TESTTHRESH)
#	$(PY) gelscore_driver.py $(SGAGG_DATA)/alignment/pab1-orthologs-sample.fa --query-gene H.sapiens --degap --min-size $(TESTSIZE) --threshold $(TESTTHRESH) 
#	$(PY) gelscore_driver.py $(DATA)/taf15.fa --min-size $(TESTSIZE) --threshold $(TESTTHRESH) --out $(DATA)/taf15-scores.txt
#	$(PY) gelscore_driver.py $(GENOMES)/scer_orf_trans.fasta --query-gene NPL3 --min-size $(TESTSIZE) --threshold $(TESTTHRESH) --out $(DATA)/npl3-scores.txt
#	$(PY) gelscore_driver.py $(GENOMES)/scer_orf_trans.fasta --query-gene PAB1 --min-size $(TESTSIZE) --threshold $(TESTTHRESH) --out $(DATA)/pab1-scores.txt

col:
	$(PY) assign-region-values.py $(DATA)/scer-regions.txt $(DATA)/scer-region-p0.01-eigen.txt "$(DATA)/entropy-distributions/entropy-*mers.txt" --adjust --colorspace rgb --out $(DATA)/scer-region-cols2.txt
#	$(PY) assign-region-values.py $(DATA)/scer-region-subset.txt $(DATA)/scer-region-p0.01-eigen.txt "$(DATA)/entropy-distributions/entropy-*mers.txt" --adjust --colorspace rgb --out $(DATA)/scer-region-subset-cols.txt

draw:
#	$(PY) draw-regions.py $(GENOMES)/scer_orf_trans.fasta $(DATA)/scer-region-subset-cols.txt --query-gene PAB1 --query-gene PBP1 --query-gene SUP35 --query-gene RNQ1 --query-gene NEW1 --query-gene URE2 --query-gene NPL3 --query-gene YRA1 --query-gene NAB2 --threshold 1 --out $(FIG)/scer-subset-domains.svg
	$(PY) draw-regions.py $(DATA)/ddx4.fa $(DATA)/scer-region-subset-cols.txt --query-gene PAB1 --query-gene PBP1 --query-gene SUP35 --query-gene RNQ1 --query-gene NEW1 --query-gene URE2 --query-gene NPL3 --query-gene YRA1 --query-gene NAB2 --threshold 1 --out $(FIG)/scer-subset-domains.svg

coltest:
	$(PY) assign-region-values.py $(DATA)/scer-regions-test.txt $(DATA)/scer-region-p0.01-eigen.txt "$(DATA)/entropy-distributions/entropy-*mers.txt" --colorspace ycbcr --out $(DATA)/scer-cols-test.txt

drawtest:
	$(PY) draw-regions.py $(GENOMES)/scer_orf_trans.fasta $(DATA)/scer-region-cols2.txt  \
		--query-gene PAB1  --query-gene PBP1 --query-gene NPL3 --query-gene PUB1 --query-gene NGR1 \
		--query-gene NRD1 --query-gene NAB3 --query-gene PUF3 --query-gene MSB1 --query-gene CAF120 \
		--query-gene SUP35 --query-gene RNQ1 --query-gene URE2 --query-gene NEW1 --query-gene SWI1 \
		--query-gene SIS1 --query-gene YDJ1 --query-gene APJ1 \
		--query-gene CDC19 --query-gene PGK1 --query-gene TDH3 --query-gene FBA1 \
		--query-gene SSA1 --query-gene SSA2 --query-gene SSA3 --query-gene SSA4 \
		--query-gene CPR6 --query-gene SSE1 --query-gene HSC82 --query-gene HSP82 \
		--query-gene GUS1 --query-gene MES1 --query-gene ARC1 \
		--query-gene SRO9 --query-gene DED1 --query-gene YRA1 \
		--query-gene HSF1 \
		--threshold 0.01 --out $(FIG)/scer-domains.svg
#		--query-gene IST2 --query-gene RPG1 --query-gene SRP40 --query-gene FUN30 \

pab: segpab colpab drawpab

segpab:
#	$(PY) gelscore_driver.py $(SGAGG_DATA)/alignment/pab1-orthologs-sample.fa --min-size 15 --degap --threshold 0.05 --out $(DATA)/pab1-orthologs-regions.txt
	$(PY) gelscore_driver.py $(SGAGG_DATA)/alignment/pab1-domains/pab1-RRRNC.fa --min-size 15 --degap --threshold 0.05 --out $(DATA)/pab1-RRRNC-regions.txt

colpab:
	$(PY) assign-region-values.py $(DATA)/pab1-orthologs-regions.txt $(DATA)/scer-region-p0.01-eigen.txt "$(DATA)/entropy-distributions/entropy-*mers.txt" --adjust --out $(DATA)/pab-region-pvals.txt

drawpab:
	$(PY) draw-regions.py $(SGAGG_DATA)/alignment/pab1-orthologs-sample.fa  $(DATA)/pab-region-pvals.txt --degap --label-width 5.0 --threshold 1 --out $(FIG)/pab-domains.svg
	
ddx: segddx colddx drawddx

segddx:
	$(PY) gelscore_driver.py $(DATA)/ddx4.fa --min-size 15 --degap --threshold 0.01 --out $(DATA)/ddx4-regions.txt

colddx:
	$(PY) assign-region-values.py $(DATA)/ddx4-regions.txt $(DATA)/scer-region-p0.01-eigen.txt "$(DATA)/entropy-distributions/entropy-*mers.txt" --adjust --out $(DATA)/ddx4-region-pvals.txt

drawddx:
	$(PY) draw-regions.py $(DATA)/ddx4.fa  $(DATA)/ddx4-region-pvals.txt --degap --label-width 5.0 --threshold 0.05 --out $(FIG)/ddx4-domains.svg
	

gen:
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 5 -n 100000 --upper-window-size 25 --out $(DATA)/entropy-distributions/entropy &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 26 -n 100000 --upper-window-size 50 --out $(DATA)/entropy-distributions/entropy &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 51 -n 100000 --upper-window-size 75 --out $(DATA)/entropy-distributions/entropy &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 76 -n 100000 --upper-window-size 100 --out $(DATA)/entropy-distributions/entropy &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 101 -n 100000 --upper-window-size 125 --out $(DATA)/entropy-distributions/entropy &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 126 -n 100000 --upper-window-size 150 --out $(DATA)/entropy-distributions/entropy &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 151 -n 100000 --upper-window-size 175 --out $(DATA)/entropy-distributions/entropy &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 176 -n 100000 --upper-window-size 200 --out $(DATA)/entropy-distributions/entropy &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 201 -n 100000 --upper-window-size 225 --out $(DATA)/entropy-distributions/entropy &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 226 -n 100000 --upper-window-size 250 --out $(DATA)/entropy-distributions/entropy &
	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 75 -n 100000 --out $(DATA)/entropy-distributions/entropy-75mers.txt &
#	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 298 -n 100000 --upper-window-size 300 --out $(DATA)/entropy-distributions/entropy &

gentest:
	$(PY) generate-window-entropy.py $(GENOMES)/scer_orf_trans.fasta 10 -n 100 --upper-window-size 15 --out $(DATA)/test

clean:
	rm -rf $(DATA)/test-*mers.txt

seg:
	$(DEV)/seg/seg.exe $(DATA)/fus.fa 20 2.3 2.5 > $(DATA)/seg-fus.txt
	$(DEV)/seg/seg.exe $(DATA)/taf15.fa 20 2.3 2.5 > $(DATA)/seg-taf15.txt

