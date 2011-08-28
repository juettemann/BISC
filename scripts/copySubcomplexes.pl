#!/usr/bin/perl

use warnings;
use strict;
use lib '/home/tmp/BISC/libs/';
use Data::Dumper;
use Pod::Usage;
use Pdbaa;
use Sql;
use Parallel::ForkManager;

pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

use constant BISCHOM => '/home/tmp/BISC/update/output/BiscHom/BiscHom.interacting'; 
use constant BISCHET => '/home/tmp/BISC/update/output/BiscHet/BiscHet.interacting'; 

use constant PDB  => '/home/tmp/BISC/update/subcomplexes/pdb/'; 
use constant PISA => '/home/tmp/BISC/update/subcomplexes/pisa/'; 

use constant PISA_MODIFIED => '/home/tmp/BISC/update/pisa/structures/modified/';
use constant PDB_STRUCTURES => '/home/tmp/BISC/update/pdb/structures/';


#use constant T_PDB_HOM  => '/home/tmp/BISC/update/upload/BiscHom/subcomplexes/pdb/';
#use constant T_PISA_HOM => '/home/tmp/BISC/update/upload/BiscHom/subcomplexes/pisa/';
#use constant T_PDB_HET  => '/home/tmp/BISC/update/upload/BiscHet/subcomplexes/pdb/';
#use constant T_PISA_HET => '/home/tmp/BISC/update/upload/BiscHet/subcomplexes/pisa/';
use constant T_PDB_HOM   => '/home/tmp/BISC/update/upload/structures/BiscHom/pdb/';
use constant T_PISA_HOM  => '/home/tmp/BISC/update/upload/structures/BiscHom/pisa/';
use constant T_PDB_HET   => '/home/tmp/BISC/update/upload/structures/BiscHet/pdb/';
use constant T_PISA_HET  => '/home/tmp/BISC/update/upload/structures/BiscHet/pisa/';

use constant T_COORD => '/home/tmp/BISC/update/upload/coordinateFiles/';

our $pdbaa = new Pdbaa;
our $sql = new Sql;
our $pm = new Parallel::ForkManager(4);

&main();
sub main(){

   copySubcomplexes(BISCHOM, T_PISA_HOM, T_PDB_HOM);
   copySubcomplexes(BISCHET, T_PISA_HET, T_PDB_HET);
}

sub copySubcomplexes(){
   my ($Lfile,$pisaTarget,$pdbTarget) = @_;
   my ($dbh,$Lsource,$Ltarget,$complex);

   my $HRia = $pdbaa->parseInteractions($Lfile);

   foreach my $id (keys %{$HRia->{dir}}){
#$pm->start and next; 
#1fnt:A:B
      my @id = split(/:/,$id);
      my $pdbID  = $id[0];
      my $chain1 = $id[1];
      my $chain2 = $id[2];

      my $subDir = substr($pdbID,1,2) . '/';
      my $complex = $pdbID . '_' . $chain1 . $chain2 . '.pdb';

      my $LsourcePisa = PISA . $subDir . $complex;
      my $LsourcePdb  = PDB  . $subDir . $complex;

      my $LtargetPisa = $pisaTarget . $complex;
      my $LtargetPdb  = $pdbTarget  . $complex;

      my $pattern = $pdbID . '.pisa.modified.*';
      my $ARfiles = $pdbaa->parseDir(PISA_MODIFIED,  $pattern);

      if(scalar @{$ARfiles} > 1) {warn PISA_MODIFIED . $pattern; print Dumper($ARfiles); }
#For some reason the binary file has vanished
# Check if chains are available in PBD & PISA, then redo
      
      if( (not -T $LsourcePisa) && (not -T $LsourcePdb) ){
         warn "$LsourcePisa and $LsourcePdb\nare missing. Shoot me.\n";
         my $Lpdb  = PDB_STRUCTURES . $pdbID . '.pdb';
         my $Lpisa = PISA_MODIFIED . $ARfiles->[0];
         my ($HRpdbFile)  = $pdbaa->parsePDBatom($Lpdb);#print Dumper($HRpdbFile); 
         my ($HRpisaFile) = $pdbaa->parsePDBatom($Lpisa );#print Dumper($HRpdbFile); 

         if( (exists $HRpdbFile->{0}->{$chain1}) && (exists $HRpdbFile->{0}->{$chain2})){
            printComplex($LsourcePdb, $HRpdbFile, $chain1, $chain2);           
         }
         if( (exists $HRpisaFile->{0}->{$chain1}) && (exists $HRpisaFile->{0}->{$chain2})){
            printComplex($LsourcePisa, $HRpisaFile, $chain1, $chain2);           
         }

      }
#copies the modified structure (1fnt.pisa.modified.1.pdb
      if($ARfiles->[0]){
         my $LsourceCoord = PISA_MODIFIED . $ARfiles->[0];
         my $LtargetCoord = T_COORD . $ARfiles->[0];
            my $args = "cp $LsourceCoord $LtargetCoord";
         if ((-T $LsourceCoord) && (not -T $LtargetCoord) ){
            system($args) == 0 or warn "$args";
         }
         undef($LsourceCoord);
         undef($LtargetCoord);
      }
      print "$LsourcePisa\n"; 
      if ((-T $LsourcePisa) && (not -T $LtargetPisa) ){
         my $args = "cp $LsourcePisa $LtargetPisa";
        system($args) == 0 or warn "$args";
      }
      undef($LsourcePisa);
      undef($LtargetPisa);
      if ( (-T $LsourcePdb) && (not -T $LtargetPdb) ){
         my $args = "cp $LsourcePdb $LtargetPdb";
         system($args) == 0 or warn "$args";
      }
#$pm->finish;
   }
#$pm->wait_all_children;
}
sub Het(){
   my ($dbh,$Lsource,$Ltarget,$complex);

   my $HRia = $pdbaa->parseInteractions(BISCHET);

   foreach my $id (keys %{$HRia->{dir}}){
      my @id = split(/:/,$id);
      my $pdbID = $id[0];
      my $chain1 = $id[1];
      my $chain2 = $id[2];

      my $subDir = substr($pdbID,1,2) . '/';
      my $complex = $pdbID . '_' . $chain1 . $chain2 . '.pdb';

      my $LsourcePisa = PISA . $subDir . $complex;
      my $LsourcePdb  = PDB  . $subDir . $complex;
      my $LtargetPdb = T_PDB_HET .  $complex;
      my $LtargetPisa = T_PISA_HET .  $complex;

      if( (not -T $LsourcePisa) && (not -T$LsourcePdb) ){warn "Pisa: $LsourcePisa\nPdb: $LsourcePdb";}

      my $pattern = $pdbID . '.pisa.modified.*';
      my $ARfiles = $pdbaa->parseDir(PISA_MODIFIED,  $pattern);
      if(scalar @{$ARfiles} != 1) {warn PISA_MODIFIED.$pattern; print Dumper($ARfiles); }
#For some reason the binary file has vanished
# Check if chains are available in PBD & PISA, then redo
      if( (not -T $LsourcePisa) && (not -T $LsourcePdb) ){
         print "Missing:\n$LsourcePisa\n$LsourcePdb\n$chain1,$chain2\n\n"; 

         my $Lpdb = PDB_STRUCTURES . $pdbID . '.pdb';
         my $Lpisa = PISA_MODIFIED . $ARfiles->[0];
         print "$Lpdb\n$Lpisa\n\n"; 
         my ($HRpdbFile) = $pdbaa->parsePDBatom(PISA_MODIFIED . $ARfiles->[0]);#print Dumper($HRpdbFile); 
         my ($HRpisaFile) = $pdbaa->parsePDBatom($Lpdb );#print Dumper($HRpdbFile); 

         if( (exists $HRpdbFile->{0}->{$chain1}) && (exists $HRpdbFile->{0}->{$chain2})){

            printComplex($LsourcePdb, $HRpdbFile, $chain1, $chain2);           
         }
         if( (exists $HRpisaFile->{0}->{$chain1}) && (exists $HRpisaFile->{0}->{$chain2})){
            printComplex($LsourcePdb, $HRpisaFile, $chain1, $chain2);           
         }

      }
      if(scalar @{$ARfiles} > 1) {warn $pattern; print Dumper($ARfiles); die; }
      if($ARfiles->[0]){
         my $LsourceCoord = PISA_MODIFIED . $ARfiles->[0];
         my $LtargetCoord = T_COORD . $ARfiles->[0];
         if ((-T $LsourceCoord) && (not -T $LtargetCoord) ){
            my $args = "cp $LsourceCoord $LtargetCoord";
            system($args) == 0 or warn "$args";
         }
         undef($LsourceCoord);
         undef($LtargetCoord);
      }
#if(not -T $LsourcePisa){die "$LsourcePisa"}
      if ((-T $LsourcePisa) && (not -T $LtargetPisa) ){
         my $args = "cp $LsourcePisa $LtargetPisa";
         system($args) == 0 or die "$args";
      }
      undef($LsourcePisa);
      undef($LtargetPisa);
      if ( (-T $LsourcePdb) && (not -T $LtargetPdb) ){
         my $args = "cp $LsourcePdb $LtargetPdb";
         system($args) == 0 or die "$args";
      }
   }
}
# print ATOM section for complex file
sub printComplex($\%$){
   my $complex = shift @_;
   my $HRpdbFile = shift @_;
   my $chain1 = shift @_;
   my $chain2 = shift @_;

   open(FH,'>',$complex) or die "Can not open/access '$complex'\n$!";
      print FH $HRpdbFile->{'0'}->{$chain1};
      print FH "TER\n"; 
      print FH $HRpdbFile->{'0'}->{$chain2};
      print FH "END";
   close(FH);
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

=begin comment
      my $subDir = substr($pdbID,1,2) . '/';
         $complex = $pdbID . '_' . $chain1 . $chain2 . '.pdb';
      if(($pdbInterface > 200) && ($pisaInterface > 200)){
         $Lsource = PISA . $subDir . $complex;
         $Ltarget = T_PISA_HET;
      }
      elsif($pisaInterface == -1){
         $Lsource = PDB . $subDir  . $complex;
         $Ltarget = T_PDB_HET;
      }
      elsif($pdbInterface == -1){
         $Lsource = PISA . $subDir . $complex;
         $Ltarget = T_PISA_HET;
      }
      if(not -e $Ltarget.$complex){
         my $args = "cp $Lsource $Ltarget";
         system($args) == 0 or die "$args";
      }
   }
   $sql->disconnectBisc(\$dbh); 
