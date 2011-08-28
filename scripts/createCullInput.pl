#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pdbaa;
use Sql;

our $pdbaa = new Pdbaa;
our $sql = new Sql;
our $stamp = $pdbaa->stampDay();

&main();
#
#  Read in new interactons, megre with old ones (from MySQL db, sot highest chains top, print lut, cull
#
#
use constant PDBAA      => '/home/tmp/BISC/update/pdbaa.new';
use constant PRE_CULL_HOM_PDB => '/home/tmp/BISC/update/output/BiscHom/PdbPreCullBiscHom';
use constant PRECULLHETPDB => '/home/tmp/BISC/update/output/BiscHet/PdbPreCullBiscHet';

use constant PRE_CULL_HOM_PISA => '/home/tmp/BISC/update/output/BiscHom/PisaPreCullBiscHom';
use constant PRECULLHETPISA => '/home/tmp/BISC/update/output/BiscHet/PisaPreCullBiscHet';

use constant CULLINHOM  => '/home/tmp/BISC/update/output/BiscHom/piscesInputBiscHom';
use constant CULLINHET  => '/home/tmp/BISC/update/output/BiscHet/piscesInputBiscHet';

use constant REFCHAINS => '/home/tmp/BISC/scripts/pdb2pisaChains.log'; 


sub hom(){
   my ($HRpdbaa,$HRsize,$HRrefChains) = @_;
   my $HRpreCullBiscHomPdb  = $pdbaa->parseInteractions(PRE_CULL_HOM_PDB);
   my $HRpreCullBiscHomPisa = $pdbaa->parseInteractions(PRE_CULL_HOM_PISA);
   my $HRsizeHom = &merge($HRpreCullBiscHomPdb,$HRpreCullBiscHomPisa,$HRsize);#print Dumper($HRsizeHom); 
   
   &cullHom($HRsizeHom,$HRpdbaa,$HRrefChains);
}
sub het(){
   my ($HRpdbaa,$HRsize,$HRrefChains) = @_;
   my $HRpreCullBiscHetPdb  = $pdbaa->parseInteractions(PRECULLHETPDB);#print Dumper($HRpreCullBiscHetPdb); die;
   my $HRpreCullBiscHetPisa = $pdbaa->parseInteractions(PRECULLHETPISA);
   my $HRsizeHet = &merge($HRpreCullBiscHetPdb,$HRpreCullBiscHetPisa,$HRsize); #print Dumper($HRsizeHet); die;
   &cullHet($HRsizeHet,$HRpdbaa,$HRrefChains);

}

sub main(){
#my $LbiscHom  = '/home/s0571283/BISC/output/preCullBiscHom';
   my ($dbh);
   $sql->connectBisc(\$dbh);
      my $HRsize = &numberofChains($dbh);
   $sql->disconnectBisc(\$dbh);

#   my $HRold = &getOldBisc($dbh);# print Dumper($HRold);die;  
   my $HRpdbaa = $pdbaa->parseDunbrackFasta(PDBAA);#print Dumper($HRdunbrack->{'sequence'});die;
   my $HRrefChains = &parsePisa2PdbMapping(REFCHAINS);
   &hom($HRpdbaa,$HRsize,$HRrefChains);
   &het($HRpdbaa,$HRsize,$HRrefChains);

}

#Join PDB and PISA IDs
sub merge(\%\%\%){
   my ($HRpdb,$HRpisa,$HRsize) = @_;
   my %result;

   foreach my $id (keys %{$HRpdb->{pdbID}}){
      my $pdbID = $HRpdb->{pdbID}->{$id};
      my $size = $HRsize->{$pdbID};
      $size = 0 if(! $size);
      $result{'size'}{$pdbID} = $size;
      $result{'ids'}{$pdbID}{$id}++;
   }
   
   foreach my $id (keys %{$HRpisa->{pdbID}}){
      my $pdbID = $HRpisa->{pdbID}->{$id};
      my $size = $HRsize->{$pdbID};
      $size = 0 if(! $size);
      $result{'size'}{$pdbID} = $size;
      $result{'ids'}{$pdbID}{$id}++;
   }
   return(\%result);

}


# get the number of each chain

sub numberofChains(){
   if(scalar @_ != 1){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 1 needed\n";die;}
   my ($dbh) =  @_;
   my %result;
   my $sth = $dbh->prepare(q{
         select pdbID, chainsPDB, proteinPisa from pdbPisaChains;
         }) or die $dbh->errstr;
   $sth->execute() or die  $dbh->errstr;
   while(my ($pdbID,$chainsPdb,$proteinPisa)  = $sth->fetchrow_array()){
         my $pdbIDlc = lc($pdbID);
         if(not exists $result{$pdbIDlc}){
            if (!$proteinPisa){$proteinPisa = 0}
            ($proteinPisa >= $chainsPdb) ? ($result{$pdbIDlc} = $proteinPisa) :  ($result{$pdbIDlc} = $chainsPdb);
         }
   }
   return(\%result);
}


sub cullHet(\%\%){
   if(scalar @_ != 3){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 3 needed\n";die;}
   my ($HRsizeHet,$HRpdbaa,$HRrefChains) = @_;
   my (@cullin,%written);
      foreach my $pdbIDlc (reverse sort {$HRsizeHet->{'size'}->{$a} <=> $HRsizeHet->{'size'}->{$b}} keys %{$HRsizeHet->{'size'}} ){
         foreach my $id (sort keys %{$HRsizeHet->{'ids'}->{$pdbIDlc}} ){
                  
               $id =~ /(\w{4}):(\w):(\w)/;
            my $pdbIDlc = $1;
            my $chain1 = $2;
            my $chain2 = $3;
            if (not defined $pdbIDlc){print "ID P: $id\n"; }
            if (not defined $chain1){print "ID C1: $id\t$chain1\n"; }
            if (not defined $chain2){print "ID C2:$id\t$chain2\n"; }
            my $pdbIDuc = uc($pdbIDlc);
            my $seq1 = $HRpdbaa->{'fasta'}->{$pdbIDuc.$chain1};
            my $seq2 = $HRpdbaa->{'fasta'}->{$pdbIDuc.$chain2};
            if(!$seq1){
               my $refChain1 = $HRrefChains->{$pdbIDlc}->{$chain1};
#      print "1: $pdbID\t$refChain1\n"; die;
               $seq1 = $HRpdbaa->{'sequence'}->{$pdbIDuc.$refChain1};
            }
            if(!$seq2){
               my $refChain2 = $HRrefChains->{$pdbIDlc}->{$chain2};
#      print "2: $pdbID\t$refChain2\t$chain2\n";
#     print Dumper($HRrefChains->{$pdbIDlc});  die;
               $seq2 = $HRpdbaa->{'sequence'}->{$pdbIDuc.$refChain2};
            }
            if( not exists $written{$pdbIDlc}{$seq1}){ push(@cullin,$seq1)}
            $written{$seq1}++;
            if( not exists $written{$pdbIDlc}{$seq2}){ push(@cullin,$seq2)}
            $written{$seq2}++;
         }
      }
      my $out = join("\n",@cullin);
      open(FH,'>',CULLINHET) or die "Can not open/access '" . CULLINHET .  "'\n$!"; 
         print FH $out; 
      close(FH);
      my $HRseq = $pdbaa->parseFastaOrdered(CULLINHET);
      my (%ref,@seq);
      foreach my $c (sort {$HRseq->{$a} <=> $HRseq->{$b}} keys %{$HRseq}){
         my $seq = $HRseq->{$c}->{sequence};
         next if( exists $ref{$seq});
         $ref{$seq}++;
         my $fasta = $HRseq->{$c}->{fasta};
         push(@seq,$fasta);
      }
      open(FH,'>',CULLINHET) or die "Can not open/access '" . CULLINHET .  "'\n$!"; 
         print  FH join("\n", @seq), "\n";
      close(FH);

}
#
sub cullHom(\%\%){
   if(scalar @_ != 3){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 3 needed\n";die;}
   my ($HRsizeHom,$HRpdbaa,$HRrefChains) = @_;
   my (%written,@cullin);
      foreach my $pdbIDlc (reverse sort {$HRsizeHom->{'size'}->{$a} <=> $HRsizeHom->{'size'}->{$b}} keys %{$HRsizeHom->{'size'}} ){
         if (not defined $HRsizeHom->{'size'}->{$pdbIDlc}){ next}; 
         foreach my $id (sort keys %{$HRsizeHom->{'ids'}->{$pdbIDlc}} ){

            $id =~ /(\w{4}):(\w):(\w)/;
            my $pdbIDlc = $1;
#            if($pdbIDlc eq '3e2j'){die} 
            next if(exists $written{$pdbIDlc});
            my $chain1 = $2;
            my $chain2 = $3;
            my $pdbIDuc = uc($pdbIDlc);
            my $seq1 = $HRpdbaa->{'fasta'}->{$pdbIDuc.$chain1};
            if(!$seq1){
               my $refChain1 = $HRrefChains->{$pdbIDlc}->{$chain1};
#      print "1: $pdbID\t$refChain1\n"; die;
               $seq1 = $HRpdbaa->{'sequence'}->{$pdbIDuc.$refChain1};
            }
            if( not exists $written{$seq1}){ push(@cullin,$seq1)}
            $written{$seq1}++;
         }
      }
      my $out = join("\n",@cullin);
      open(FH,'>',CULLINHOM) or die "Can not open/access '" . CULLINHOM .  "'\n$!"; 
         print FH  $out; 
      close(FH);
#quick fix
      my $HRseq = $pdbaa->parseFastaOrdered(CULLINHOM);
      my (%ref,@seq);

      foreach my $c (sort {$HRseq->{$a} <=> $HRseq->{$b}} keys %{$HRseq}){
         my $seq = $HRseq->{$c}->{sequence};
         next if( exists $ref{$seq});
         $ref{$seq}++;
         my $fasta = $HRseq->{$c}->{fasta};
         push(@seq,$fasta);
      }
      open(FH,'>',CULLINHOM) or die "Can not open/access '" . CULLINHOM .  "'\n$!"; 
         print  FH join("\n", @seq), "\n";
      close(FH);

}


=begin comment
#PDB ID | PDB reference chain | PISA original chain | PISA new chain
10gs  A  A  A
=cut
sub parsePisa2PdbMapping($){
   my %result;
   my ($Lfile) = @_;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^#/);
         chomp($line);
         my @line = split("\t",$line); 
         my $pdbID = $line[0];
         my $ref = $line[1];
         my $org = $line[2];
         my $new = $line[3];
         $result{$pdbID}{$new} = $ref;

      }
   close(FH);
return(\%result);

}

=begin comment
#pdbID:chain1:chain2 sequenceIdentity interfaceSize asaSingleChain1 asaSingleChain2 asaComplexChain1 asaComplexChain2
2BDG:A:B 100.0 434   9995  9972  9531  9567
2PGL:A:B 100.0 487   13557 13591 13075 13098
2R9H:A:B 100.0 3571  19376 19219 15791 15661
=cut

=begin comment
sub getOldBisc($){
   if(scalar @_ != 1){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 1 needed\n";die;}
   my ($dbh) =  @_;
   my (%old);
   my $sth = $dbh->prepare(q{
      select pdbID, chainID1, chainID2 from biscHet;
      }) or die $dbh->errstr;
   $sth->execute() or die  $dbh->errstr;

   while(my ($pdbID,$chain1,$chain2)  = $sth->fetchrow_array()){
      $old{$pdbID}{$chain1}++;
      $old{$pdbID}{$chain2}++;
   }
   
   $sth = $dbh->prepare(q{
      select pdbID, chainID1, chainID2 from biscHom;
      }) or die $dbh->errstr;
   $sth->execute() or die  $dbh->errstr;

   while(my ($pdbID,$chain1,$chain2)  = $sth->fetchrow_array()){
      $old{$pdbID}{$chain1}++;
      $old{$pdbID}{$chain2}++;
   }

   return(\%old);
}
sub pdbaaNumberChainsPDB(\%){
   if(scalar @_ != 1){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 1 needed\n";die;}
   my ($HRpdbaa) =  @_;
   my %result;
   foreach my $pdbID (keys %{$HRpdbaa->{'list'}}){
      $result{$pdbID} = scalar(@{$HRpdbaa->{'list'}->{$pdbID}});
   }
   return(\%result);
}
sub createCullInputBiscHet(\%\%\%$){
   if(scalar @_ != 4){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 4 needed\n";die;}
   my ($HRpdbaa, $HRpreCullBiscHet, ) =  @_;
   my $stamp = 
   my $sthpdb = $dbh->prepare(q{
         select pdbID, chainID1, chainID2 from biscHet;
         }) or die $dbh->errstr;
   $sthpdb->execute() or die  $dbh->errstr;
   while ( my($pdbID,$chain1,$chain2) = $sthpdb->fetchrow_array() ) {
      $HRpreCullBiscHet->{'line'}->{$pdbID}{$chain1}{$chain2} = 'old interaction';
   }

   $sthpdb->finish();

   open(FH,'>',CULLINHET . $stamp) or die "Can not open/access '" . CULLINHET . $stamp . "'\n$!"; 
#complexes with most chains top
   foreach my $pdbID (reverse sort {$HRsize->{$a} <=> $HRsize->{$b}  }keys %{$HRsize}){
      
      if(exists $HRpreCullBiscHet->{'line'}->{$pdbID}){
         foreach my $chain1 (sort keys %{$HRpreCullBiscHet->{'line'}->{$pdbID}}){
            if(not exists $HRpdbaa->{'fasta'}->{$pdbID.$chain1}){ print "$pdbID\t$chain1\n";next;}
            print FH $HRpdbaa->{'fasta'}->{$pdbID.$chain1}."\n"; 
            foreach my $chain2 (sort keys %{$HRpreCullBiscHet->{'line'}->{$pdbID}->{$chain1}}){
               if(not exists $HRpdbaa->{'fasta'}->{$pdbID.$chain2}){ print "$pdbID\t$chain2\n";next;}
               print FH $HRpdbaa->{'fasta'}->{$pdbID.$chain2}."\n"; 
            }
         }
      }
   }
}

sub createCullInputBiscHom(\%\%){
   if(scalar @_ != 4){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 4 needed\n";die;}
   my ($HRpdbaa, $HRpreCullBiscHom, $HRsize, $dbh) = @_;

   my $sthpdb = $dbh->prepare(q{
         select pdbID, chainID1, chainID2 from biscHet;
         }) or die $dbh->errstr;
   $sthpdb->execute() or die  $dbh->errstr;
   while ( my($pdbID,$chain1,$chain2) = $sthpdb->fetchrow_array() ) {
      $HRpreCullBiscHom->{'line'}->{$pdbID}{$chain1}{$chain2} = 'old interaction';
   }

   $sthpdb->finish();

   open(FH,'>',CULLINHOM . $stamp) or die "Can not open/access '".CULLINHOM . $stamp ."'\n$!"; 
#complexes with most chains top
   foreach my $pdbID (reverse sort {$HRsize->{$a} <=> $HRsize->{$b}  }keys %{$HRsize}){
      
      if(exists $HRpreCullBiscHom->{'line'}->{$pdbID}){
         foreach my $chain1 (sort keys %{$HRpreCullBiscHom->{'line'}->{$pdbID}}){
            print FH $HRpdbaa->{'fasta'}->{$pdbID.$chain1}."\n"; 
            last;
         }
      }
   }
}
=head1 NAME

createCullInput.pl


=head1 SYNOPSIS

createCullInput.pl

=head1 AUTHOR

Thomas Juettemann <juettemann@gmail.com>


=head1 ARGUMENTS

=over 8

=item none

=back


=head1 OPTIONS

=over 8

=item none

=back

=head1 DESCRIPTION

B<createCullInput.pl> creates input file for PISCES culling 

=head1 BUGS

=over 4

=item A few (10) pdbIDs are not found
 
=back

=head1 TODO

=over 4

=item Nothing

=back

=head1 UPDATES

=cut


