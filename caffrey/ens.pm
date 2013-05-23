#!/usr/bin/perl
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use strict;
my $registry = 'Bio::EnsEMBL::Registry';

sub connect{
my $registry = 'Bio::EnsEMBL::Registry'; $registry->load_registry_from_db(
-host => 'ensembldb.ensembl.org',
-user => 'anonymous'
);}

sub getSpliceAdaptorSpecies{
  my $species = shift;
  my $sa = $registry->get_adaptor($species,"core","slice");
  return $sa;
}

sub printGeneCanonInfo(){
my ($gene,$speciesName) = @_;
    my $trans = $gene->canonical_transcript();
    my $chr = $gene->slice->seq_region_name();
    my $id = $trans->stable_id();
    my $length = $trans->length();
    my $can3UTR = $trans->three_prime_utr();
    my $can5UTR = $trans->five_prime_utr();
    my $threeLen = 0; my $fiveLen  = 0;
    if (defined $can3UTR){$threeLen = $can3UTR->length();}
    if (defined $can5UTR){$fiveLen = $can5UTR->length();}
    if ($threeLen > 0 or $fiveLen > 0){
     print join(",",$id,$speciesName,$chr,$length,$threeLen,$fiveLen)."\n";
    }
}

sub pOutUtrInfo{
my ($sa,$speciesName) = @_;
my @genes = map { @{$_->get_all_Genes}  } @{$sa->fetch_all('chromosome')};
map { &printGeneCanonInfo($_,$speciesName) } @genes;
}


1;
