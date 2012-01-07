#!/usr/bin/perl

use 5.10.0;

use strict;
use warnings;

use Log::Log4perl qw(get_logger);
use Data::Dumper;
use XML::LibXML;
use Carp;

use feature qw(say);

use PdbMeta;

#use feature qw(switch say);

my $file = 'pdb1fin.ent';
my $meta  = PdbMeta->new();


&main();

sub main(){
   my $parser = XML::LibXML->new();
   my $cfg    = $parser->parse_file('config.xml');
   my $entDir = $cfg->findvalue('/cfg/path/pdb_ent_dir');

  my $xray = identifyCrystalStructures($cfg->findvalue('/cfg/path/pdbaa'));
  say 'Read '. keys(%{$xray}).' xray structures';

  foreach my $pdbID (sort keys %{$xray}){
    my $entFile = buildPdbFilePath($pdbID, 'ent');
       $entFile = $entDir . $entFile;

    my $amountStructures = $meta->getNumberOfBiomolecules($entFile);
    confess "Error\n$entFile: $amountStructures" if(!$amountStructures);
# test how many chains are in if only one
# if more than one, test if chain names are duplicated, rewrite
# copy results to new folder
    next if(scalar(@{$amountStructures}) == 1);
    say $entFile;
    my $remark350 = $meta->parseRemark350($entFile);
    print Dumper($remark350); die;


# check if assemblies with different amount of chains exist.
# get the one with more chains if
# if several with the same amount of chains exists, use the one with lowest energy
# 
# read chain section, check if chain IDs have been used multiple times
# if so, rename keep relation
# split, identify interacting files
# homo/hetero distinction
# cluster

#print Dumper($amountStructures); 
#print Dumper($remark350); 
  }
sub _testForSingleChain {
#use pdbaa 
}

}
###############################################################################
#                                   _findLargestAssembly
###############################################################################
# Identify the largest assembly with the lowest free energy in PDB ent files
# In many cases, the file contains this remark:
# REMARK 300 BIOMOLECULE: 1, 2
# The numbers seem to concure to the PDB files available in the biounit sub-
# directory. Parsing REMARK 350, that should be available in all files,
# the assembly/structure with the most chains/lowest free energy is identified
# It is also checked if author and software agree on the biounit

###############################################################################
sub _findLargestAssembly {
  my ($remark350) = @_;
  my $last;
  $last->{chains} = 0;
  $last->{energy} = 1000;
  $last->{biounit} = 0;
# Check if both author and software determined values exists for molecules
# Compare if they are equal
# Compare if they are a) more chains than former molecule b) delta G is lower
  foreach my $unit (sort keys %{$remark350}){

    if( (exists $remark350->{$unit}->{author}) && 
        (exists $remark350->{$unit}->{software}) ){

      if($remark350->{$unit}->{author} != $remark350->{$unit}->{software}){
        die "Different, check it" 
      }

      if($remark350->{$unit}->{software} > $last->{chains}){
        $last->{chains}  = $remark350->{$unit}->{software};
        $last->{energy}  = $remark350->{$unit}->{energy};
        $last->{biounit} = $remark350->{$unit}->{biounit};
      }
      if($remark350->{$unit}->{software} = $last->{chains}){
        if(defined $remark350->{$unit}->{energy} and 
            $remark350->{$unit}->{energy} < $last->{energy}){

          $last->{energy}  = $remark350->{$unit}->{energy};
          $last->{biounit} = $remark350->{$unit}->{biounit};
        }
      }
      next;
    }
# If only software exists, only replace if more chains or equal chains and lower gibbs
# ToDo: PISA/PQS
    if(exists $remark350->{$unit}->{software}){
      if($remark350->{$unit}->{software} > $last->{chains}){
        $last->{chains} = $remark350->{$unit}->{software};
        $last->{energy}   = $remark350->{$unit}->{energy};
        $last->{biounit}  = $remark350->{$unit}->{biounit};
      }
      if($remark350->{$unit}->{software} =  $last->{chains}){
        if(defined $remark350->{$unit}->{energy} and 
            $remark350->{$unit}->{energy} < $last->{energy}){
          $last->{chains} = $remark350->{$unit}->{software};
          $last->{energy}   = $remark350->{$unit}->{energy};
          $last->{biounit}  = $remark350->{$unit}->{biounit};
        }
      }
      next;
    }
#
    if(exists $remark350->{$unit}->{author}){
      if($remark350->{$unit}->{author} > $last->{chains}){
        $last->{chains}   = $remark350->{$unit}->{author};
        $last->{energy}   = $remark350->{$unit}->{energy};
        $last->{biounit}  = $remark350->{$unit}->{biounit};
      }
      next;
    }
  }
}

# 1fnt-> fn/1nft
sub buildPdbFilePath {
  my ($pdbID,$type) = @_;
  $type = 'normal' if (!$type);

  my $subDir = substr($pdbID,1,2);
  my $path;
  $path = "$subDir/pdb$pdbID.ent" if ($type eq 'ent');
  $path = "$subDir/$pdbID.pdb"    if ($type eq 'normal');
  return($path);

}
# read pdbaa, keep crystal structures
#>200LA 164 XRAY
#PDB+chain ID given, return only PDB ID in lowercase
sub identifyCrystalStructures {
  my ($file) = @_;
  say $file;
  my $result;
  open(my $fh,'<', $file) or confess "Can not open/access '$file'\n$!";
    while(my $line = <$fh>){
      next unless($line =~ />(\w{5})\s\d*\sXRAY/);
      my $pdbID = substr($1,0,4);
      $pdbID = lc($pdbID);
      $result->{$pdbID}++;
    }
  close($fh);
return($result);


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


