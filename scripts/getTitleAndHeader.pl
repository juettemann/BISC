#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pdbaa;
use Sql;

our $sql = new Sql;
our $pdbaa = new Pdbaa;

use constant INTERACTINGHOM => '/home/tmp/BISC/update/output/BiscHom/BiscHom.interacting';
use constant INTERACTINGHET => '/home/tmp/BISC/update/output/BiscHet/BiscHet.interacting';
use constant DIRECTORIES    => '/home/tmp/BISC/update/pdb/structures/';
use constant OUTPUT         => '/home/tmp/BISC/update/output/Bisc/BiscToTitleHeader';

&main();

sub main(){
#my $LinteractingFile = shift @ARGV;
   my %result;
   my $HRinteractingFile1    = $pdbaa->parseInteractions(INTERACTINGHOM);
   my $HRinteractingFile2    = $pdbaa->parseInteractions(INTERACTINGHET);
   my %merged     = (%{$HRinteractingFile1->{'ids'}},%{$HRinteractingFile2->{'ids'}});
#   my %merged     = %{$HRinteractingFile1->{'ids'}};


   open(LOG,'>',$0.'.log') or die "Can not open/access '".$0.".log'.\n$!\n";
   select((select(LOG), $|=1)[0]); 
   print LOG "Total PDB: ".keys (%merged)."\n";

   foreach my $pdbID (keys %merged){
      print LOG "Started: $pdbID\n"; 
      my $pdbFile = DIRECTORIES . $pdbID .'.pdb';


      if(-e $pdbFile){
         my $header = $pdbaa->getPdbHeader($pdbFile);
         my $title = $pdbaa->getPdbTitle($pdbFile);
         if(not defined $header){$header = 'Not available';}
         if(not defined $title){$title = 'Not available';}
         $result{$pdbID}{'header'} = $header;
         $result{$pdbID}{'title'} = $title;
      }
      else{warn "$pdbFile not found"}
      print LOG "Finished: $pdbID\n"; 
   } 
   close(LOG);
   open(FH,'>',OUTPUT) or die "Can not open/access '".OUTPUT."'\n$!"; 
   foreach my $pdbID (keys %result){
      print FH "$pdbID\t$result{$pdbID}{'title'}\t$result{$pdbID}{'header'}\n";
   } 
   close(FH);
   my $dbh;
   $sql->connectBisc(\$dbh);

   $dbh->do('truncate pdb') or die $dbh->errstr;

   my $sth = $dbh->prepare(q{
         insert into pdb (pdbID,title,header) values (?,?,?);
         }) or die $dbh->errstr;

   foreach my $pdbID (sort keys %result){
      my $title = $result{$pdbID}{'title'};
      my $header = $result{$pdbID}{'header'};

      $sth->execute($pdbID, $title,$header) or die  $dbh->errstr;
   }

   $sql->disconnectBisc(\$dbh); 

   
}
