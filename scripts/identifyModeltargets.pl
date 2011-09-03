#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pdbaa;
use Sql;

our $pdbaa = new Pdbaa;
our $Sql   = new Sql;

use constant PDBAA             => '/home/tmp/BISC/update/pdbaa.23Aug2010';

use constant INTERACTINGHOM    => '/home/tmp/BISC/update/output/BiscHom/BiscHom.interacting';
use constant INTERACTINGHET    => '/home/tmp/BISC/update/output/BiscHet/BiscHet.interacting';

use constant INTERACTINGINTACT      => '/home/tmp/BISC/intact/intact.buck.ppi';
use constant INTERACTINGHPRD        => '/home/tmp/BISC/hprd/hprd.buck.ppi';
use constant INTERACTINGBIOGRID     => '/home/tmp/BISC/biogrid/biogrid.buck.ppi';
use constant INTERACTINGPLASMODIUM  => '/home/tmp/BISC/input/plasmodium/interaction_data_gene_info_051109.csv';


use constant BLASTOUTINTACTHOM  => '/home/tmp/BISC/update/output/BiscHom/buck.intact.Aug28.biscHom.Aug23.blast.out';
use constant BLASTOUTHPRDHOM    => '/home/tmp/BISC/update/output/BiscHom/buck.hprd.release9.biscHom.Aug23.blast.out';
use constant BLASTOUTBIOGRIDHOM => '/home/tmp/BISC/update/output/BiscHom/buck.biogrid.release3-0-67.biscHom.Aug23.blast.out';

use constant BLASTOUTINTACTHET  => '/home/tmp/BISC/update/output/BiscHet/buck.intact.Aug28.biscHet.Aug23.blast.out';
use constant BLASTOUTHPRDHET    => '/home/tmp/BISC/update/output/BiscHet/buck.hprd.release9.biscHet.Aug23.blast.out';
use constant BLASTOUTBIOGRIDHET => '/home/tmp/BISC/update/output/BiscHet/buck.biogrid.release3-0-67.biscHet.Aug23.blast.out';

use constant BLASTOUTHOMPFAL  => '/home/tmp/BISC/blast/blast-2.2.20/plasmodium20090811VsBiscHom27Mar2009.blast.out';
use constant BLASTOUTHETPFAL  => '/home/tmp/BISC/blast/blast-2.2.20/plasmodium20090811VsBiscHet27Mar2009.blast.out';

use constant OUTHPRDHET       => '/home/tmp/BISC/update/output/BiscHet/buck.none.modeltargets23August2010.hprd.biscHet';
use constant OUTBIOGRIDHET    => '/home/tmp/BISC/update/output/BiscHet/buck.none.modeltargets23August2010.biogrid.biscHet';
use constant OUTINTACTTHET    => '/home/tmp/BISC/update/output/BiscHet/buck.none.modeltargets23August2010.intact.biscHet';
use constant OUTPLASMODIUMHET => '/home/tmp/BISC/update/output/BiscHet/buck.modeltargets23August2010.plasmodium.biscHet.noEvalLimit';

use constant OUTHPRDHOM       => '/home/tmp/BISC/update/output/BiscHom/buck.none.modeltargets23August2010.hprd.biscHom';
use constant OUTBIOGRIDHOM    => '/home/tmp/BISC/update/output/BiscHom/buck.none.modeltargets23August2010.biogrid.biscHom';
use constant OUTINTACTTHOM    => '/home/tmp/BISC/update/output/BiscHom/buck.none.modeltargets23August2010.intact.biscHom';
use constant OUTPLASMODIUMHOM => '/home/tmp/BISC/update/output/BiscHom/modeltargets23August2010.plasmodium.biscHom.noEvalLimit';


#use constant INTERACTINGINTACT      => '/home/s0571283/BISC/html/viva/intact.viva.ppi';
#use constant OUTINTACTTHOM    => '/home/s0571283/BISC/html/viva/modeltargets.viva.intact.biscHom';
#use constant OUTINTACTTHET    => '/home/s0571283/BISC/html/viva/modeltargets.viva.intact.biscHet';

############## Change table names (modeltargetsHet) ##################################

&main();
sub main(){
   
   my $HRbiscHetInteractions = $pdbaa->parseInteractions(INTERACTINGHET) ;   
   my $HRbiscHomInteractions = $pdbaa->parseInteractions(INTERACTINGHOM) ;   
   my $HRpdbaa               = $pdbaa->parseDunbrackFasta(PDBAA);
   
   &hprd($HRbiscHetInteractions,$HRbiscHomInteractions,$HRpdbaa);
   &biogrid($HRbiscHetInteractions,$HRbiscHomInteractions,$HRpdbaa);
   &intact($HRbiscHetInteractions,$HRbiscHomInteractions,$HRpdbaa);
#   &plasmodium($HRbiscHetInteractions,$HRbiscHomInteractions,$HRpdbaa);
}
sub intact(\%\%\%){
   my ($HRbiscHetInteractions,$HRbiscHomInteractions,$HRpdbaa) = @_;
   my $table = 'modeltargetsIntact';
   my $HRintactInteractions  = $pdbaa->parseIntActinteractions(INTERACTINGINTACT) ;   #print Dumper($HRintactInteractions); die;
#   print "Start BLAST parse\n"; sleep(5);

   my $HRintactBlastOutHet   = $pdbaa->parseUniprotLikeVsBisc(BLASTOUTINTACTHET) ;
   my $HRtargetsHet = &identifyTargets($HRintactInteractions,$HRbiscHetInteractions->{'line'},$HRintactBlastOutHet,$HRpdbaa);#print Dumper($HRtargetsHom);
   &printResults($HRtargetsHet,OUTINTACTTHET);

   my $HRintactBlastOutHom   = $pdbaa->parseUniprotLikeVsBisc(BLASTOUTINTACTHOM) ;
   my $HRtargetsHom = &identifyTargets($HRintactInteractions,$HRbiscHomInteractions->{'line'},$HRintactBlastOutHom,$HRpdbaa);#print Dumper($HRtargetsHom);
   &printResults($HRtargetsHom,OUTINTACTTHOM);

#   &populateModeltargets($HRtargetsHet,$table);

}
sub hprd(\%\%\%){
   my ($HRbiscHetInteractions,$HRbiscHomInteractions,$HRpdbaa) = @_;
   my $table = 'modeltargetsHPRD';
#works for all, rename
   my $HRhprdInteractions  = $pdbaa->parseIntActinteractions(INTERACTINGHPRD) ;   #print Dumper($HRhprdInteractions); die;

  my $HRhprdBlastOutHet   = $pdbaa->parseUniprotLikeVsBisc(BLASTOUTHPRDHET) ;
   my $HRtargetsHet = &identifyTargets($HRhprdInteractions,$HRbiscHetInteractions->{'line'},$HRhprdBlastOutHet,$HRpdbaa);#print Dumper($HRtargetsHom);
   &printResults($HRtargetsHet,OUTHPRDHET);

   my $HRhprdBlastOutHom   = $pdbaa->parseUniprotLikeVsBisc(BLASTOUTHPRDHOM) ;
   my $HRtargetsHom = &identifyTargets($HRhprdInteractions,$HRbiscHomInteractions->{'line'},$HRhprdBlastOutHom,$HRpdbaa);#print Dumper($HRtargetsHom);
   &printResults($HRtargetsHom,OUTHPRDHOM);

#   &populateModeltargets($HRtargetsHet,$table);

}

sub biogrid(\%\%\%){
   my ($HRbiscHetInteractions,$HRbiscHomInteractions,$HRpdbaa) = @_;
   my $table = 'modeltargetsBIOGRID';
   my $HRbiogridInteractions  = $pdbaa->parseIntActinteractions(INTERACTINGBIOGRID) ;   #print Dumper($HRbiogridInteractions); die;

   my $HRbiogridBlastOutHet   = $pdbaa->parseUniprotLikeVsBisc(BLASTOUTBIOGRIDHET) ;
   my $HRtargetsHet = &identifyTargets($HRbiogridInteractions,$HRbiscHetInteractions->{'line'},$HRbiogridBlastOutHet,$HRpdbaa);#print Dumper($HRtargetsHom);
   &printResults($HRtargetsHet,OUTBIOGRIDHET);

   my $HRbiogridBlastOutHom   = $pdbaa->parseUniprotLikeVsBisc(BLASTOUTBIOGRIDHOM) ;
   my $HRtargetsHom = &identifyTargets($HRbiogridInteractions,$HRbiscHomInteractions->{'line'},$HRbiogridBlastOutHom,$HRpdbaa);#print Dumper($HRtargetsHom);
   &printResults($HRtargetsHom,OUTBIOGRIDHOM);

#   &populateModeltargets($HRtargetsHet,$table);

}

sub plasmodium(\%\%\%){
   if(scalar @_ != 3){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 3 needed\n";die;}
   my ($HRbiscHetInteractions,$HRbiscHomInteractions,$HRpdbaa) = @_;
   my $table = 'modeltargetsPlasmodium';
   my $HRplasmodiumInteractions  = $pdbaa->parsePlasmosiumInteractions(INTERACTINGPLASMODIUM) ;   #print Dumper($HRbiogridInteractions); die;
   my $HRplasmodiumBlastOutHom   = $pdbaa->parsePlasmodiumVsBisc(BLASTOUTHOMPFAL) ;
   my $HRplasmodiumBlastOutHet   = $pdbaa->parsePlasmodiumVsBisc(BLASTOUTHETPFAL) ;

   my $HRtargetsHet = &identifyTargets($HRplasmodiumInteractions,$HRbiscHetInteractions->{'line'},$HRplasmodiumBlastOutHet,$HRpdbaa);#print Dumper($HRtargetsHom);
   &printResults($HRtargetsHet,OUTPLASMODIUMHET);
   my $HRtargetsHom = &identifyTargets($HRplasmodiumInteractions,$HRbiscHomInteractions->{'line'},$HRplasmodiumBlastOutHom,$HRpdbaa);#print Dumper($HRtargetsHom);
   &printResults($HRtargetsHom,OUTPLASMODIUMHOM);
}
#"bait_orf","product_name","prey_orf","product_name"
#"MAL13P1.120","splicing factor, putative","MAL13P1.202","conserved Plasmodium protein, unknown function"
#"MAL13P1.120","splicing factor, putative","MAL8P1.111","JmjC domain containing protein"
sub parsePlasmosiumInteractions($){
   if(scalar @_ != 1){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 1 needed\n";die;}
   my ($Lfile) =  @_; 
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^"bait/);
         chomp($line);
         $line =~ s/^"//;
         $line =~ s/"$//;
         my @line = split(/","/,$line);
         $result{$line[0]}{$line[2]}{'bait'} = $line[1];
         $result{$line[0]}{$line[2]}{'prey'} = $line[3];
      }   
   close(FH);
   return(\%result);
}

#for debuging
sub printBlast(){
   my ($HRblastResult) =  @_;
   foreach my $upid (sort {$HRblastResult->{$a} <=> $HRblastResult->{$b}} keys %{$HRblastResult}){
      foreach my $c (sort  keys %{$HRblastResult->{$upid}}){
         foreach my $pdbAndChainID (sort  keys %{$HRblastResult->{$upid}->{$c}}){
            if (($pdbAndChainID eq '2NRUB') && ($upid eq 'Q6H6T1')){
               print "$upid\n";
               print "$c\n";
               print "$pdbAndChainID\n"; 
               print Dumper $HRblastResult->{$upid}->{$c}->{$pdbAndChainID}; 
            }
         }
      }
   }
}

sub identifyTargets(\%\%\%\%){
   if(scalar @_ != 4){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 4 needed\n";die;}
   die "Passed a non hashref as argument 1" if (ref($_[0]) ne 'HASH');
   die "Passed a non hashref as argument 2" if (ref($_[1]) ne 'HASH');
   die "Passed a non hashref as argument 3" if (ref($_[2]) ne 'HASH');
   die "Passed a non hashref as argument 4" if (ref($_[3]) ne 'HASH');


   my $FGinteractions = shift @_;#foreach my $key (keys %{$FGinteractions}){foreach my $key2 (keys %{$FGinteractions->{$key}}){print "$key\t$key2\n";}} die;#print Dumper($FGinteractions); die;
   my $HRbiscInteractions   = shift @_;
   my $HRBlastresults       = shift @_;#print Dumper($HRBlastresults);die; 
   my $HRpdbaa              = shift @_;
#print Dumper($HRBlastresults); die;
   my (%result,%temp,%missingSeq);
   my $counter = 0;
#iterate through all functional genomics interactions
LOOP:
   foreach my $id1 (keys %{$FGinteractions}){
      foreach my $id2 (keys %{$FGinteractions->{$id1}}){       #if(length($id2)<4){print "ID1:$id1\nId2$id2\n\n"; }
#print "-1\n"; 
# see if interaction matched in BLAST
         if( (exists $HRBlastresults->{$id1}) && (exists $HRBlastresults->{$id2})){ 
#print "0\n"; 
# go top->bottom through blast results for each interaction partner
            foreach my $rank1 (    sort {$HRBlastresults->{$id1}->{$a} <=> $HRBlastresults->{$id1}->{$b} }  keys %{$HRBlastresults->{$id1}}){
               foreach my $rank2 ( sort {$HRBlastresults->{$id2}->{$a} <=> $HRBlastresults->{$id2}->{$b} }  keys %{$HRBlastresults->{$id2}}){
#get PDB and chain IDs for each interaction partner
                  foreach my $pdbAndChainID1 (   sort keys %{ $HRBlastresults->{$id1}->{$rank1}}){
#       print "$HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{db}\n";die; 
                     next if($HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{db} !~ /HUMAN$/);
                     foreach my $pdbAndChainID2 (sort keys %{ $HRBlastresults->{$id2}->{$rank2}}){
                        next if($HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{db} !~ /HUMAN$/);
                        my $pdbAndChainID1lc = lc(substr($pdbAndChainID1,0,4)) . substr($pdbAndChainID1,4,1);
                        my $pdbAndChainID2lc = lc(substr($pdbAndChainID2,0,4)) . substr($pdbAndChainID2,4,1);

                        my $pdbAndChainID1uc = uc(substr($pdbAndChainID1,0,4)) . substr($pdbAndChainID1,4,1);
                        my $pdbAndChainID2uc = uc(substr($pdbAndChainID2,0,4)) . substr($pdbAndChainID2,4,1);
#consistency check     
                        if(not defined $HRpdbaa->{'sequence'}->{$pdbAndChainID1uc})
                        {print "Missing sequence: $pdbAndChainID1lc\n"; die;
                           $missingSeq{$pdbAndChainID1uc}++;
                        }
                        if(not defined $HRpdbaa->{'sequence'}->{$pdbAndChainID2uc}){
                           print "Missing sequence: $pdbAndChainID2uc\n";
                           $missingSeq{$pdbAndChainID2uc}++;
                        }

                        my $pdbID1   = substr($pdbAndChainID1lc,0,4);
                        my $pdbID2   = substr($pdbAndChainID2lc,0,4);
                        
                        unless($pdbID1 eq $pdbID2){next;}
#print "1\n";                         
                        my $length1  = length($HRpdbaa->{'sequence'}->{$pdbAndChainID1uc});
                        my $length2  = length($HRpdbaa->{'sequence'}->{$pdbAndChainID2uc});
                        my $chain1 = substr($pdbAndChainID1lc,4,1);
                        my $chain2 = substr($pdbAndChainID2lc,4,1);
                        
                        if($chain1 eq $chain2){next;}
#print "2\n";                         
                        
                        unless( (exists $HRbiscInteractions->{$pdbID1}->{$chain1}->{$chain2}) ||
                                (exists $HRbiscInteractions->{$pdbID1}->{$chain2}->{$chain1})){next}
#print "3\n";                         
                        unless ($HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'aliLength'} > 50){next;}  
#print "4\n";                         
                        unless ($HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'aliLength'} > 50){next;}   
#                        unless ($HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'%id'} > 25){next;}  
#                        unless ($HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'%id'} > 25){next;}  
#print "5\n";                         

                        unless ($length1 >= ($HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'slength'}/2)){next;}
#print "6\n";                         
                        unless ($length2 >= ($HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'slength'}/2)){next;}
#print "7\n";                         

                        my $eval1 = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'evalue'}; 
                        my $eval2 = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'evalue'}; 
                        
                        if( (exists $temp{$id1}{$id2}) || (exists $temp{$id2}{$id1}) ){next}
                        else{$temp{$id1}{$id2}++}
#print "8\n";                         

                        my @chains = ($chain1,$chain2);
                        @chains = sort(@chains);
                        if(($chain1 eq $chains[0]) && ($chain2 eq $chains[1])){

                           $result{$counter}{'pdbID'}      = $pdbID1;
                           $result{$counter}{'chain1'}     = $chain1;
                           $result{$counter}{'chain2'}     = $chain2;
                           $result{$counter}{'uniprotID1'} = $id1;
                           $result{$counter}{'uniprotID2'} = $id2;
                           $result{$counter}{'evalue1'}    = $eval1;
                           $result{$counter}{'evalue2'}    = $eval2;
                           $result{$counter}{'%id1'}       = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'%id'};
                           $result{$counter}{'%id2'}       = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'%id'};
                           $result{$counter}{'aliLength1'} = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'aliLength'};
                           $result{$counter}{'aliLength2'} = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'aliLength'};
                           $result{$counter}{'slength1'}   = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'slength'};
                           $result{$counter}{'slength2'}   = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'slength'};
                           $result{$counter}{'q.start1'}   = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'q.start'};
                           $result{$counter}{'q.start2'}   = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'q.start'};
                           $result{$counter}{'q.end1'}     = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'q.end'};
                           $result{$counter}{'q.end2'}     = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'q.end'};
                           $result{$counter}{'s.start1'}   = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'s.start'};
                           $result{$counter}{'s.start2'}   = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'s.start'};
                           $result{$counter}{'s.end1'}     = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'s.end'};
                           $result{$counter}{'s.end2'}     = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'s.end'};
                           $counter++;
#                           next LOOP;
                        }
                        elsif(($chain1 eq $chains[1]) && ($chain2 eq $chains[0]) ){
                           $result{$counter}{'pdbID'}      = $pdbID1;
                           $result{$counter}{'chain1'}     = $chain2;
                           $result{$counter}{'chain2'}     = $chain1;
                           $result{$counter}{'uniprotID1'} = $id2;
                           $result{$counter}{'uniprotID2'} = $id1;
                           $result{$counter}{'evalue1'}    = $eval2;
                           $result{$counter}{'evalue2'}    = $eval1;
                           $result{$counter}{'%id1'}       = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'%id'};
                           $result{$counter}{'%id2'}       = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'%id'};
                           $result{$counter}{'aliLength1'} = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'aliLength'};
                           $result{$counter}{'aliLength2'} = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'aliLength'};
                           $result{$counter}{'slength1'}   = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'slength'};
                           $result{$counter}{'slength2'}   = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'slength'};
                           $result{$counter}{'q.start1'}   = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'q.start'};
                           $result{$counter}{'q.start2'}   = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'q.start'};
                           $result{$counter}{'q.end1'}     = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'q.end'};
                           $result{$counter}{'q.end2'}     = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'q.end'};
                           $result{$counter}{'s.start1'}   = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'s.start'};
                           $result{$counter}{'s.start2'}   = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'s.start'};
                           $result{$counter}{'s.end1'}     = $HRBlastresults->{$id2}->{$rank2}->{$pdbAndChainID2}->{'s.end'};
                           $result{$counter}{'s.end2'}     = $HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}->{'s.end'};
                           $counter++;
#                           next LOOP;

                        }
                        else{
                           print "$pdbID1\t$chain1\t$chain2\n"; 
                           print Dumper(\@chains);  
                           print Dumper($HRBlastresults->{$id1}->{$rank1}->{$pdbAndChainID1}) ."\n"; 
                           print Dumper($HRBlastresults->{$id2}->{$rank1}->{$pdbAndChainID2}) ."\n"; 
                           die;
                        }
#print "\n\n\n";                         

                     }
                  }
               }
            }  
         }
      }
   }
#   print Dumper(\%missingSeq); 
#print Dumper(\%result);die; 
   return(\%result);
}
sub printResults($){
   my ($HRmodeltargets,$out) =  @_;

   print "$out\n"; 

   open(FH,'>',$out ) or die "Can not open/access '".$out."'\n$!\n"; 
      print FH "#pdbID, chain1, chain2, uniprotID1, uniprotID2, evalue1, evalue2, \%identity1, \%identity2,aliLength1,aliLength2,slength1,slength2,";
      print FH "q.start1,q.start2,q.end1,q.end2,s.start1,s.start2,s.end1,s.end2\n"; 
      foreach my $counter (keys %{$HRmodeltargets}){
         print FH $HRmodeltargets->{$counter}->{'pdbID'}	     . "\t";
         print FH $HRmodeltargets->{$counter}->{'chain1'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'chain2'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'uniprotID1'} . "\t";
         print FH $HRmodeltargets->{$counter}->{'uniprotID2'} . "\t";
         print FH $HRmodeltargets->{$counter}->{'evalue1'}    . "\t"; 
         print FH $HRmodeltargets->{$counter}->{'evalue2'}    . "\t"; 
         print FH $HRmodeltargets->{$counter}->{'%id1'}	     . "\t";
         print FH $HRmodeltargets->{$counter}->{'%id2'}       . "\t"; 
         print FH $HRmodeltargets->{$counter}->{'aliLength1'} . "\t";
         print FH $HRmodeltargets->{$counter}->{'aliLength2'} . "\t";
         print FH $HRmodeltargets->{$counter}->{'slength1'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'slength2'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'q.start1'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'q.start2'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'q.end1'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'q.end2'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'s.start1'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'s.start2'}	  . "\t";
         print FH $HRmodeltargets->{$counter}->{'s.end1'}     . "\t";
         print FH $HRmodeltargets->{$counter}->{'s.end2'}     . "\n";

      } 
   close(FH);
}


sub parseHprdInteractions($){
   my $Lfile = shift @_;
   my %result;
   print "$Lfile\n"; 

   open(FH,'<',$Lfile) or die "Can not open/access '" . $Lfile ."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^#/);
         chomp($line);
         my @line = split(/\t/,$line);#print Dumper(\@line);die; 
            $result{$line[2]}{$line[3]}++;
      }
   close(FH);
   return(\%result);
}
sub parseBiogridInteractions($){
   my ($Lfile) =  @_;
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '" . $Lfile ."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^#/);
         chomp($line);
         my @line = split(/\t/,$line);
            $result{$line[4]}{$line[5]}++;
      }
   close(FH);
   return(\%result);
}
#####Same method in Pdbaa#########
#CHEBI:15367 uniparc:UPI00000738EA 
#CHEBI:15414 uniprotkb:O43159 uniprotkb:P17095-1 uniprotkb:P17095-2 uniprotkb:P17096-1 uniprotkb:Q96LA8 uniprotkb:Q9987
#OUTPUT . IntActInteractions
sub parseInteractions($){
   my $Lfile = shift @_;
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '" . $Lfile ."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^#/);
         chomp($line);
         my @line = split(/\t/,$line);
            $result{$line[2]}{$line[3]}++;
      }
   close(FH);
   return(\%result);
}
