#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pdbaa;


our $pdbaa = new Pdbaa;

&main();

sub main(){
   my $Lin = '/home/tmp/BISC/update/output/BiscHet/piscesInputBiscHet';
   my $HRseq = $pdbaa->parseFastaOrdered($Lin);
   my (%ref,@seq);

   foreach my $c (sort {$HRseq->{$a} <=> $HRseq->{$b}} keys %{$HRseq}){
      my $seq = $HRseq->{$c}->{sequence};
      next if( exists $ref{$seq});
      $ref{$seq}++;
      my $fasta = $HRseq->{$c}->{fasta};
      push(@seq,$fasta);
   }
      print  join("\n", @seq), "\n";

}
