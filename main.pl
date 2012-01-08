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
use PdbParser;

#use feature qw(switch say);

my $file = 'pdb1fin.ent';
my $meta  = PdbMeta->new();
my $parser = PdbParser->new();

&main();

sub main(){
   my $xml_parser = XML::LibXML->new();
   my $cfg    = $xml_parser->parse_file('config.xml');
   my $entDir = $cfg->findvalue('/cfg/path/pdb_ent_dir');

  my $xray = identifyCrystalStructures($cfg->findvalue('/cfg/path/pdbaa'));
  say 'Read '. keys(%{$xray}).' xray structures';

  foreach my $pdbID (sort keys %{$xray}){
next unless($pdbID eq '2biw');
    my $entFile = buildPdbFilePath($pdbID, 'ent');
       $entFile = $entDir . $entFile;
#say $entFile;
    #How many biomolecules in the file?
    my $remark350 = $meta->parseRemark350($entFile);
    if(keys $remark350 > 1){
      my ($size, $bioUnitID, $softwareInfo) =
        _getBioMolecule($remark350, $entFile);
      say "Size: $size, Unit: $bioUnitID, $softwareInfo";
#How many chains in the biomolecule with the most chains
      if($size > 1){
        _buildBioUnitPath($cfg, $pdbID, $unit);
        die;
      }
    }
#    print Dumper($remark350); die;
# Deprecated: my $amountStructures = $meta->getNumberOfBiomolecules($entFile);
#    confess "Error\n$entFile: $amountStructures" if(!$amountStructures);
# test how many chains are in if only one
# if more than one, test if chain names are duplicated, rewrite
# copy results to new folder

#   if(scalar(@{$amountStructures}) == 1) {
#    my($chains, $chains2array_pos) = $parser->readChains($entFile);
#    if (scalar(@{$chains}) > 1){
# check if chain IDs are duplicated, replace
# write new file to new directory

#    }
#print Dumper($chains);
#    print Dumper($chains2array_pos);
#    die;
#   }
  }
}
###############################################################################
#
###############################################################################
# /home/tmp/BISC/pdb/data/biounit/coordinates/divided/bi/2biw.pdb1.gz
sub _buildBioUnitPath {
  my ($cfg, $pdbID, $unit) = @_;
  my $bioUnitDir = $cfg->findvalue('/cfg/path/pdb_biounit_dir');

  my $subDir = substr($pdbID,1,2);
  my $path = "$subDir/$pdbID.pdb.$unit";
  my $gz_path = "$subDir/$pdbID.pdb.$unit.gz";
  if (-e $gz_path){
    my $args = "gunzip $gz_path";
    system($args) == 0 or die $args;
  }
  return($path);

}
###############################################################################
#                             _FindBestStructure
###############################################################################
# The largest biomolecule will be chosen
sub _getBioMolecule {
  my ($remark350, $entFile) = @_;

  my $tmp;

  foreach my $molNo (sort keys %{$remark350}){

    my ($software, $softwareInfo) = _buildSoftwareInfo($remark350->{$molNo});
    my ($author);

    if(exists $remark350->{$molNo}->{author}){
      $author = $remark350->{$molNo}->{author};
    }
# In software we trust
    if(defined $software && defined $author ){
# check if author and software agree
      if ($author != $software){
        warn "$entFile difference between author and software";
        if($author>$software){
          push( @{$tmp->{$author}},"$molNo:a");
        }
        else {
          push( @{$tmp->{$software}},"$molNo:$softwareInfo");
        }
      }
      else {
        push( @{$tmp->{$software}},"$molNo:$softwareInfo");
      }
    }
    elsif(exists $remark350->{$molNo}->{software}){
      push( @{$tmp->{$software}},"$molNo:$softwareInfo");
    }
    elsif(exists $remark350->{$molNo}->{author}){
      push( @{$tmp->{$author}},"$molNo:a");
    }
    else{
      die "PDB always has a surprise left...";
    }
  }
  my ($size, $biomolecule, $softwareInfo) = _chooseBestBioMolecule($tmp);
}

###############################################################################
#
###############################################################################
sub _chooseBestBioMolecule {
  my ($tmp) = @_;
# Find largest structure, pick on with this hierachy:
# 1st choice: PISA, 2nd: Author, 3rd: PQS
  foreach my $size (reverse sort {$a <=> $b} keys %{$tmp}){
#say "Size: $size";
#print Dumper($tmp->{$size});
    for my $bioMol (@{$tmp->{$size}}){
      if($bioMol =~ /(\d+):s:(PISA:.*)$/){
        return ($size, $1, $2);
      }
    }
#No PISA found, try author
    for my $bioMol (@{$tmp->{$size}}){
      if($bioMol =~ /(\d+):a/){
        return ($size, $1, 0);
      }
    }
#No author found, try PQS
    for my $bioMol (@{$tmp->{$size}}){
#      say "PQS: $bioMol";
      if($bioMol =~ /(\d+):s:(PQS:.*)$/){
        return ($size, $1, $2);
      }
    }
    die "PDB always has a surprise left...";
  }
  last;
}

###############################################################################
#
###############################################################################
# Get as much information from the used software as possible
sub _buildSoftwareInfo {
  my ($record) = @_;
  my ($software,
      $method,
      $energy,
      $buried,
      $surface,
      $softwareInfo
     ) = (0, 0, 0, 0, 0, 0);

  if(exists $record->{software}){
    $software = $record->{software};
    if (exists $record->{method}){
      $method = $record->{method}
    }
    if (exists $record->{energy}){
      $energy = $record->{energy}
    }
    if (exists $record->{buried}){
      $buried = $record->{buried}
    }
    if (exists $record->{surface}){
      $surface = $record->{surface}
    }
    $softwareInfo = "s:$method:$energy:$buried:$surface";
  }
  return($software, $softwareInfo);
}

###############################################################################
#
###############################################################################
sub _renameDuplicates {
  my ($chains2array_pos, $chains) = @_;

  # leave only chains as potential new ones that are not present yet
  my @range = ('A'..'Z', 'a'..'z', 0..9);
  my @chainsInPdb = sort keys($chains2array_pos);
  for my $chain(@chainsInPdb) {
    @range = grep { $_ ne $chain } @range;
  }
  # One never knows...
  return(0) if( scalar(@range) < 1);

  foreach my $chain (sort keys %{$chains2array_pos}){
    # same chain ID more than once in assembly
    my $size = scalar @{$chains2array_pos->{$chain}};
    if($size > 1){
      return (0) if($size > scalar(@range) );
      # leave 1st chain as it is, change all following
      for(my $i = 1; $i < $size; $i++) {
        my $index       = $chains2array_pos->[$i];
        my $chainRecord = $chains->[$index];
        my $newChain    = shift(@range);
        my @newRecord;
        for my $row(@{$chainRecord}) {
          say $row;
          $row = substr($row,21,1,$newChain);
          say $row;
        }
      }
    }
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


