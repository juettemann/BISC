#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use Pdbaa;
use Sql;

pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

use constant PDBAA  => '/home/tmp/BISC/update/pdbaa.april26';
use constant LOG_HOM  => '/home/tmp/BISC/update/output/BiscHom/log_pc95.log';
use constant LOG_HET  => '/home/tmp/BISC/update/output/BiscHet/log_pc95.log';

our $pdbaa = new Pdbaa;
our $sql = new Sql;

&main();
sub main(){
   my ($dbh);
   $sql->connectBisc(\$dbh);
   my $HRpdbaa = &getPdbXray();
#   print Dumper($HRpdbaa);die; 
#   &Hom($dbh,$HRpdbaa);
   &Het($dbh,$HRpdbaa);
   $sql->disconnectBisc(\$dbh); 
}
sub Hom(){
   my ($dbh,$HRxray) = @_;
   my $HR;
   my $HRlog = &parseLog(LOG_HOM);
#   print Dumper($HRlog);die; 
   my $HRseq = &getBiscSequences($dbh,'Hom');
#   print Dumper($HRseq); 
   my $HRref = &createRefList($HRxray,$HRlog,$HRseq);
   &populatePdbRef($dbh,'Hom',$HRref);
}
sub Het(){
   my ($dbh,$HRxray) = @_;
   my $HR;
   my $HRlog = &parseLog(LOG_HET);
#   print Dumper($HRlog);die; 
   my $HRseq = &getBiscSequences($dbh,'Het');
#   print Dumper($HRseq); 
   my $HRref = &createRefList($HRxray,$HRlog,$HRseq);
   &populatePdbRef($dbh,'Het',$HRref);
   
}
sub populatePdbRef(){
   my ($dbh,$FhomHet,$HRref) = @_;
   my ($sth);
   if($FhomHet eq 'Hom'){
      $dbh->do('truncate biscHomRef');
      $sth = $dbh->prepare(q{
            insert into  biscHomRef (pdbIDref, chainIDref, pdbIDreject, chainIDreject, seqid )  values (?,?,?,?,?) ;
            }) or die $dbh->errstr;
   }

   elsif($FhomHet eq 'Het'){
      $dbh->do('truncate biscHetRef');
      $sth = $dbh->prepare(q{
            insert into  biscHetRef (pdbIDref, chainIDref, pdbIDreject, chainIDreject, seqid )  values (?,?,?,?,?) ;
            }) or die $dbh->errstr;
   }
   else{die "$FhomHet"}

   foreach my $refID (sort keys %{$HRref}){
      $refID =~ /(\w{4})(\w)/;
      my $pdbIDref = $1;
      my $chainIDref = $2;
      foreach my $rejectID (keys %{$HRref->{$refID}}){
         $rejectID =~ /(\w{4})(\w)/;
         my $pdbIDreject = $1;
         my $chainIDreject = $2;
         my $seqId = $HRref->{$refID}->{$rejectID};
#        print Dumper($HRref->{$refID}); 
#        print "$rejectID\n$pdbIDreject\n$chainIDreject\n$seqId\n"; 
#        print "\n"; 
#        sleep(1);

         $sth->execute($pdbIDref,$chainIDref,$pdbIDreject,$chainIDreject,$seqId) or die  $dbh->errstr;

      }
   } 
$sth->finish();
}
sub createRefList(\%\%\%){
  my ($HRxray, $HRlog, $HRseq) = @_;
  my ($HR,$HRref);

  foreach my $id (reverse sort keys %{$HRseq}){
     my $seq = $HRseq->{$id};
     print "ID:$id\nSeq:$seq\n\n";sleep(1); 
     $HRref->{$seq}->{refID} = $id;
  }
  foreach my $id (keys %{$HRxray}){
     my $seq = $HRxray->{$id}->{SEQ};
#my $refID = $HRref->{$seq}->{refID};
#     if( !$seq or !$refID ){print Dumper($HRxray->{$id}); sleep(1); }
     if( exists $HRref->{$seq} ){
        my $refID = $HRref->{$seq}->{refID};
        if ($id ne $refID){
           $HRref->{$seq}->{reject}->{$id} = '100';
        }
     }
  }
  foreach my $seq (keys %{$HRref}){
      my $refID = $HRref->{$seq}->{refID};
      $HR->{$refID} = $HRref->{$seq}->{reject};
  } 
  foreach my $rejectID (keys %{$HRlog}){
      my $refID = $HRlog->{$rejectID}->{pdbRef};
      $HR->{$refID}->{$rejectID} = $HRlog->{$rejectID}->{seqId};
  }
return($HR);
}
#reject 3HMJC 2UV8A 99
sub parseLog(){
   my ($Llog) = @_;
   my $HR;
   open(FH,'<',$Llog) or die "Can not open/access '".$Llog."'\n$!";
      while(my $line = <FH>){
         chomp($line);
         my @line = split(/\s/,$line);
         my $reflc = lc(substr($line[1],0,4) ) . substr($line[1],4,1);
         my $rejectlc = lc(substr($line[2],0,4) ) . substr($line[2],4,1);
         my $seqId  = $line[3];
         $HR->{$rejectlc}->{pdbRef} = $reflc;
         $HR->{$rejectlc}->{seqId} = $seqId;
      }
   close(FH);
   return($HR);
}
###################################################################
# Retreive X-ray PDB entried
###################################################################
sub getPdbXray(){
   my $HR;
   my $HRpdbaa = $pdbaa->parseDunbrackFastaNew(PDBAA);
   foreach my $id (keys %{$HRpdbaa}){
         if($HRpdbaa->{$id}->{EXP} &&  ($HRpdbaa->{$id}->{EXP} eq 'XRAY') ){
            my $idlc = lc(substr($id,0,4)) . substr($id,4,1);
            $HR->{$idlc} = $HRpdbaa->{$id};

         }
   }
   return($HR);
}
###################################################################
# Retreive all sequences in BISC
###################################################################
sub getBiscSequences(){
   my ($dbh,$FhomHet) = @_;
   my ($HR,$sthBisc);
   if($FhomHet eq 'Hom'){
      $sthBisc = $dbh->prepare(q{
            select pdbID, chainID1, chainID2  from biscHom;
            }) or die $dbh->errstr;
   }
   elsif($FhomHet eq 'Het'){
      $sthBisc = $dbh->prepare(q{
            select pdbID, chainID1, chainID2  from biscHet;
            }) or die $dbh->errstr;
   }
   else{die "$FhomHet"}
      my $sthSeq = $dbh->prepare(q{
            select sequenceSeqres  from pdbSequence where pdbID = ? and chainID = ?;
            }) or die $dbh->errstr;
      $sthBisc->execute() or die  $dbh->errstr;
      while(my ($pdbID,$chain1,$chain2)  = $sthBisc->fetchrow_array()){
         $sthSeq->execute($pdbID,$chain1) or die  $dbh->errstr;
         $HR->{$pdbID.$chain1} = $sthSeq->fetchrow_array();
         $sthSeq->execute($pdbID,$chain2) or die  $dbh->errstr;
         $HR->{$pdbID.$chain2} = $sthSeq->fetchrow_array();
         
      }
   return($HR);
}
=head1 NAME

script.pl


=head1 SYNOPSIS

script.pl

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

B<script.pl> 

=head1 BUGS

=over 4

=item None known
 
=back

=head1 TODO

=over 4

=item Nothing

=back

=head1 UPDATES

=cut


