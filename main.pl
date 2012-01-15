#!/usr/bin/perl

use 5.10.0;

use strict;
use warnings;

use Log::Log4perl qw(get_logger);
use File::Path qw(make_path);
use Parallel::ForkManager;
use Data::Dumper;
use XML::LibXML;
use Carp;

use feature qw(say);

use PdbMeta;
use PdbParser;

#use feature qw(switch say);
use constant MAX_PROCESSES => 4;

my $file = 'pdb1fin.ent';
my $meta  = PdbMeta->new();
my $parser = PdbParser->new();

&main();
###############################################################################
#
###############################################################################
sub main(){
  my $pm          = new Parallel::ForkManager(MAX_PROCESSES);
  my $xml_parser  = XML::LibXML->new();
  my $cfg         = $xml_parser->parse_file('config.xml');
  my $entDir      = $cfg->findvalue('/cfg/path/pdb_ent_dir');

  my $xray = identifyCrystalStructures($cfg->findvalue('/cfg/path/pdbaa'));
  say 'Read '. keys(%{$xray}).' xray structures';

  # No REMARK 300 or 350 in 1qjb
  # 2j0u should be monomeric or dimeric, but contains 3 chains
  foreach my $pdbID (sort keys %{$xray}){
    my $pid = $pm->start and next;
    if ($pdbID =~ /^(1qjb)|(2j0u)$/) {
      $pm->finish;
      next;
    }
#next unless($pdbID eq '3q6d');
    say "Started: $pdbID";
    my $entFile = buildPdbFilePath($pdbID, 'ent');
       $entFile = $entDir . $entFile;
say $entFile;

    #How many biomolecules in the file?
    my $remark350 = $meta->parseRemark350($entFile);
    #next unless keys $remark350 > 1;
    # More than one biomolecule, we need to identify the best one
    if(keys $remark350 > 1){
      my ($size, $bioUnitID, $softwareInfo) =
        _getBioMolecule($remark350, $entFile);

      # How many chains in the biomolecule that has the most chains
      # Obviously only continue if more than 1 chain
      if($size > 1){
        _createBinaryFiles($cfg, $pdbID, $bioUnitID, $entFile);

      }
    }# more than 1 biomolecule available
    # Only 1 biomolecule
    else{
      my $bioUnitID = 1;
      _createBinaryFiles($cfg, $pdbID, $bioUnitID, $entFile);
    }
    $pm->finish;
  }
  $pm->wait_all_children;
}
###############################################################################
#                               _CreateBinaryFiles
###############################################################################
sub _createBinaryFiles {
  my ($cfg, $pdbID, $bioUnitID, $entFile) = @_;

  my $biounit = _buildBioUnitPath($cfg, $pdbID, $bioUnitID);
  my($chains, $chains2array_pos) = $parser->readChains($entFile);
  if(keys $chains2array_pos > 1){
#        print Dumper($chains);
    my $renamed = _renameDuplicates($chains2array_pos, $chains);
    die $entFile if(!$renamed);
#        print Dumper($chains2array_pos);die;
# Build path
# rename duplicates
    _splitFiles($cfg,$pdbID,$chains);
  }

}
###############################################################################
#                           RewriteChains2array_pos
###############################################################################
sub rewriteChains2array_pos {
  my ($chains2array_pos, $chains) = @_;
#emptying hash
  for (keys %{$chains2array_pos}) { delete $chains2array_pos->{$_} }

  my $size = scalar(@{$chains});
  for (my $i = 0; $i < $size ; $i++) {
    for my $line (@{$chains->[$i]}){
      my $chain = substr($line, 21, 1);
      push(@{$chains2array_pos->{$chain}}, $i);
      last;
    }
  }
}
###############################################################################
#                                 _SplitFiles
###############################################################################
sub _splitFiles {
  my ($cfg, $pdbID, $chains,$chains2array_pos) = @_;
  my $path = $cfg->findvalue('/cfg/path/splitFiles');
  my $dir = substr($pdbID,1,2);
  $path   .= "$dir/$pdbID";
  make_path($path,{verbose => 1});

  my $size = scalar(@{$chains});
#  print Dumper($chains);
  for (my $outside = 0 ; $outside < $size ; $outside++) {
    my $chainA = substr($chains->[$outside]->[0], 21, 1);
    my $file_single = "$path/$pdbID" . "_$chainA.pdb" ;
    my $a = join("\n",@{$chains->[$outside]});
    if(!-T $file_single){
      open(my $fh,'>', $file_single) or die "Can not open/access '$file_single'\n$!";
        say $fh $a;
        say $fh 'END';
      close($fh);
    }
    for (my $inside = $outside + 1 ; $inside < $size ; $inside++) {
      my $chainB = substr($chains->[$inside]->[0], 21, 1);
      my $file_complex = "$path/$pdbID" . "_$chainA"."$chainB.pdb" ;
      if(!-T $file_complex){
        my $b = join("\n",@{$chains->[$inside]});
        open(my $fh,'>', $file_complex) or die "Can not open/access '$file_complex'\n$!";
          say $fh $a;
          say $fh 'TER';
          say $fh $b;
          say $fh 'END';
        close($fh);
      }

    }
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
###############################################################################
#
###############################################################################
# /home/tmp/BISC/pdb/data/biounit/coordinates/divided/bi/2biw.pdb1.gz
###############################################################################
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
# So far the biomolecule with the most chains is chosen
# Hierachy:
# PISA, Author, PQS, unknown Software.
#
# To Do: In case of PISA, the structure with the lowest energy should be
# chosen
###############################################################################
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
    if($software && $author ){
# check if author and software agree
      if ($author != $software){
        warn "$entFile difference between author and software";
        if($author > $software){
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
  return($size, $biomolecule, $softwareInfo);
}

###############################################################################
#                           _ChooseBestBioMolecule
###############################################################################
#
#$VAR1 = {
#          '1' => [
#                   '1:s:PQS:N/A:N/A:N/A',
#                   '2:s:PQS:N/A:N/A:N/A',
#                   '3:s:PQS:N/A:N/A:N/A',
#                   '4:s:PQS:N/A:N/A:N/A'
#                 ]
#        };
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
# Unknown software used (e.g.: pdb2wdb.ent)
    for my $bioMol (@{$tmp->{$size}}){
      if($bioMol =~ /(\d+):s:/){
        return ($size, $1, 0);
      }
    }
    die "PDB always has a surprise left...";
  }
  last;
}

###############################################################################
#                             _BuildSoftwareInfo
###############################################################################
# Get as much information from the used software as possible
###############################################################################
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
#                             _RenameDuplicates
###############################################################################
# $VAR1 = {
#           'A' => [
#                    0
#                  ],
#           'B' => [
#                    1
#                  ]
#         };
###############################################################################
sub _renameDuplicates {
  my ($chains2array_pos, $chains) = @_;

  # @range contains the available new IDs
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
    next if($size == 1);
    say $chain;
    print Dumper($chains2array_pos);
    print Dumper($chains2array_pos->{$chain});
    die;

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
        die;
      }
    }
  }

  rewriteChains2array_pos($chains2array_pos, $chains);
  return(1);
}
###############################################################################
#                             BuildPdbFilePath
###############################################################################
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
###############################################################################
#                         IdentifyCrystalStructures
###############################################################################
# read pdbaa, keep crystal structures
# >200LA 164 XRAY
# PDB+chain ID given, return only PDB ID in lowercase
###############################################################################
sub identifyCrystalStructures {
  my ($file) = @_;

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


