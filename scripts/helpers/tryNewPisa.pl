package Pdbaa;


use 5.008008;
use strict;
use warnings;
use Data::Dumper;

sub parseDunbrackLog($);
sub parseDunbrackFasta($);
sub parseDunbrackFastaOld($);
sub parseWebPiscesResults($);
sub parseLocalPiscesResults($);
sub parseModeltargets();
sub parseForPdbAndChainID($);
sub createAtomSequence($$);
sub parseAstral100($);
sub parseBlast7($);
sub parseBlast7C($);
sub parseUniprotLikeVsBisc($$);
sub parseBiscVsAstralBlastOutput($);
sub parseBlastM9($);
sub parseNaccessRsa($);
sub parseNaccessAsa($);
sub parsePDBatom($);
sub parseMODELLERatom($);
sub parsePDBatomResidue($);
sub parsePDB($);
sub parsePQR($);
sub parseDirClaScoTxt($);
sub parseBiscToScop($);
sub parseDirDesScopTxtNew($);
sub parseDirDesScopTxt($);
sub parseInterfaceFile($);
sub parseInterface($);
sub parseUniprotProtText($);
sub parseFastaOrdered($);
sub parseFasta($);
sub parseUniprotFasta($);
sub parseHPRDfasta($);
sub parseSwissProtFasta($);
sub parseIntActFasta($);
sub parseHB2file($);
sub parseNcbiNamesDmp($);
sub parseMissingBiscTaxids($);
sub parseIntActToOrg($);
sub getPdbHeader($);
sub getPdbTitle($);
sub parseAstralFasta($);
sub executeNeedleWater($$$);
sub parseDir($$$);
sub stamp();
sub stampDay();
sub parseInteractions($);
sub parseBiscToOrg($);
sub parseBiscToTitleHeader($);
sub getAtomFasta;
sub parsePqrErrorFiles ($);
sub parseHprdIDMap($);
sub parseHprdPPI($);
sub parseBioGridIDs($);



sub parseDunbrackLog($){
   my ($class,$location) =  @_;
#reject 1XG0C 1XG0D 100

   my %result;

   open(LOG, "<$location") or die "Can not open/access: '$location'\n$!\n";
      while(my $line = <LOG>){
         if ($line !~ /^reject/){print STDERR "Line starts not with reject:\n$line";next;}
         chomp($line);
         $line =~ /^reject\s(\w{5})\s(\w{5})\s{1,3}(\d{1,3})$/; #external input, regex instead of split to get a error in case format changes
         my $pdbID1 = substr($1,0,4);
         my $chainID1 = substr($1,4,1);
         my $pdbID2 = substr($2,0,4);
         my $chainID2 = substr($2,4,1);
         my $seqId = $3;
         $result{$pdbID2}{$chainID2}{$pdbID1}{$chainID1}{$seqId}++;
      }
   close(LOG);

   return (\%result);
}
=begin comment
>200LA 164 XRAY 1.95 0.176 NA no LYSOZYME <UNP LYCV_BPT4> [ENTEROBACTERIA PHAGE T4]
MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAI
=cut
sub parseDunbrackFasta($){
   my ($class, $Lfile) = @_;
   my %result;
   my $reading_seq = 0;

   my ($line,$id,$sequence);
   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!\n";

LINE:
   while($line = <FH>){
      if($line =~ /^pdbid/){next;}
      chomp($line);
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1;
              my @line = split(/\s/,$line);
              my $pdbID   = substr($line,1,4);if (not defined  $pdbID){print STDERR $line;}
              my $chainID = substr($line,5,1);
              $id = $pdbID . $chainID;
              
              $result{fasta}{$id} = "$line\n";
              $result{commentLine}{$id} = $line;
              $result{pdbID}{$id} = $pdbID;
              $result{chainID}{$id} = $chainID;
              $result{ids}{$pdbID}{$chainID}++;
              push(@{$result{list}{$pdbID}},$chainID);
              $result{experiment}{$pdbID} = $line[2];        

              last SWITCH;
           };

           if ( $reading_seq && $line !~ /^>/ ) {

              $line =~ s/\s+//g;
              $result{'sequence'}{$id} .= $line;
              $sequence .= $line;
              last SWITCH;
           };

           if ( $reading_seq && $line =~ /^>/ ) {
              $sequence  =~ s/(.{80}|.{1,79}$)/$1\n/g;
              $result{'fasta'}{$id} .= $sequence;
              $sequence = '';
              $reading_seq = 0;
              redo LINE;
           };
        }

   }
   close(FH);

   return(\%result);
}

=begin comment
>200LA 164 XRAY 1.95 0.176 NA no LYSOZYME <UNP LYCV_BPT4> [ENTEROBACTERIA PHAGE T4]
MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAI
=cut
sub parseDunbrackFastaNew($){
   my ($class, $Lfile) = @_;
   my (%result,$rec,$line,$id,$seq,$fasta,$header,$pdbID,$chain,$exp,@list); 
   my $reading_seq = 0;

   open(FH,'<',$Lfile) or die "Can not open/access '$Lfile'\n$!\n";

LINE:
   while($line = <FH>){
      if($line =~ /^pdbid/){next;}
      chomp($line);
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1;
              my @line = split(/\s/,$line);
              $pdbID   = substr($line,1,4);if (!$pdbID){print STDERR $line;}
              $chain = substr($line,5,1);
              $exp = $line[2];
              $header = $line;
              $id = $pdbID . $chain;
              $fasta = "$line\n";
              last SWITCH;
           };

           if ( $reading_seq && $line !~ /^>/ ) {

              $line =~ s/\s+//g;
              $result{'seq'}{$id} .= $line;
              $seq .= $line;
              last SWITCH;
           };

           if ( $reading_seq && $line =~ /^>/ ) {
              $seq  =~ s/(.{80}|.{1,79}$)/$1\n/g;
              chomp($seq);
              $fasta .= $seq;
              $seq =~ s/\n//g;
              $rec = {
                  HEADER => $header,
                  PDBID => $pdbID,
                  CHAIN => $chain,
                  EXP => $exp,
                  FASTA => $fasta,
                  SEQ => $seq
              };
              $result{$id} = $rec;
              ($id,$header,$pdbID,$chain,$exp,$fasta,$seq) = '';

              $reading_seq = 0;
              redo LINE;
           };
        }

   }
   close(FH);

   return(\%result);


} 
# 0           1        2       3           4        5
# IDs         length   Exptl.  resolution  R-factor FreeRvalue
# 1JK9B       249      XRAY        2.900    0.22    0.26
# 1RJLC        95      XRAY        2.600    0.19    0.23

sub parseWebPiscesResults($){
   my $class = shift @_; 
   my $Lfile = shift @_;
   my %result; 

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         if($line =~ /^IDs/){next;}
         chomp($line);
         my @line    = split(/\s+/,$line);
         my $pdbID   = substr($line,0,4);
         my $chainID = substr($line,4,1);
            $result{$line[0]}{'pdbID'} = $pdbID;
            $result{$line[0]}{'chainID'} = $chainID;
            $result{$line[0]}{'length'} = $line[1];
            $result{$line[0]}{'exptl'} = $line[2];
            $result{$line[0]}{'resolution'} = $line[3];
            $result{$line[0]}{'R-factor'} = $line[4];
            $result{$line[0]}{'FreeRvalue'} = $line[5];

      }
   close(FH);
   return(\%result);
}
=cut
IDs                                              length
1DUBA                                               261
1AC6A                                               110
1IGWA                                               434
=cut

sub parseLocalPiscesResults($){
   my $class = shift @_; 
   my $Lfile = shift @_; 
   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
   while(my $line = <FH>){
      if($line =~ /^IDs/){next;}
      next if ($line =~ /^\s/);
      chomp($line);

      my @line = split(/\s+/,$line);
      my $pdbID = substr($line[0],0,4);
      my $chainID = substr($line[0],4,1);
      my $length = $line[1];

next if(!$pdbID and $length == 0);
#if(!$pdbID and $length == 0){print "$Lfile\n"}
#if(!$chainID){print "$Lfile\n"}
      $result{$pdbID}{$chainID} = $length;
   }   
   close(FH);
   return(\%result);
}

# IN (3rd line):
# 0                 1     2     3     4     5     6              7                 8     9     10    11    12    13 
# Query             Sub   %id   e.val len   bit   Org            Query             Sub   %id   e.val len   bit   Org 
# uniprotkb|Q8RTL3  2CKFD 32.53 3e-17 166   81.6  Comamonas sp.  uniprotkb|Q8RTL4  2CKFE 40.09 1e-90 439   327   Comamonas sp

# OUT:
# 2IBX' => {
#   'JK' => {
#     'organismMT2' => 'Pseudomonas pavonaceae',
#     'length1' => '55',
#     'uniprotID1' => 'uniprotkb|Q9EV84',
#     'organismMT1' => 'Pseudomonas pavonaceae',
#     'uniprotID2' => 'uniprotkb|Q9EV85',
#     'chainID2' => 'K',
#     'evalue1' => '3e-26',
#     'length2' => '61',
#     'pdbID2' => '1S0Y',
#     'evalue2' => '1e-31',
#     'pdbID1' => '1S0Y',
#     'bitScore2' => '127',
#     'chainID1' => 'J',
#     '%id2' => '100.00',
#     'bitScore1' => '109',
#     '%id1' => '100.00'
#    }
#  }


sub parseModeltargets(){
   my $class = shift @_;
   my $Lmodeltargets = shift @_;
   my %result;

   open(FH,"$Lmodeltargets") or die "Can not open/access '$Lmodeltargets'\n$!\n";
      my $c = 0;
      while(my $line = <FH>){
         chomp($line);
         my @line = split(/\t/,$line);
         my @justTheID1 = split(/\|/,$line[0]); 
         my @justTheID2 = split(/\|/,$line[7]); 
         my $pdbID1 = substr($line[1],0,4);
         my $chainID1 = substr($line[1],4,1);
         my $pdbID2 = substr($line[8],0,4);
         my $chainID2 = substr($line[8],4,1);
         if($pdbID1 ne $pdbID2){my $d = $c+1; die "Line $d PDB ID differs, the thing that should not be.\n";}
         $result{$pdbID1}{$chainID1.$chainID2}{'uniprotID1'} = $line[0];
         $result{$pdbID1}{$chainID1.$chainID2}{'justID1'} = $justTheID1[1];
         $result{$pdbID1}{$chainID1.$chainID2}{'pdbID1'} = $pdbID1;
         $result{$pdbID1}{$chainID1.$chainID2}{'chainID1'} = $chainID1;
         $result{$pdbID1}{$chainID1.$chainID2}{'%id1'} = $line[2];
         $result{$pdbID1}{$chainID1.$chainID2}{'evalue1'} = $line[3];
         $result{$pdbID1}{$chainID1.$chainID2}{'length1'} = $line[4];
         $result{$pdbID1}{$chainID1.$chainID2}{'bitScore1'} = $line[5];
         $result{$pdbID1}{$chainID1.$chainID2}{'organismMT1'} = $line[6];
         
         $result{$pdbID1}{$chainID1.$chainID2}{'uniprotID2'} = $line[7];
         $result{$pdbID1}{$chainID1.$chainID2}{'justID2'} = $justTheID2[1];
         $result{$pdbID1}{$chainID1.$chainID2}{'pdbID2'} = $pdbID2;
         $result{$pdbID1}{$chainID1.$chainID2}{'chainID2'} = $chainID2;
         $result{$pdbID1}{$chainID1.$chainID2}{'%id2'} = $line[9];
         $result{$pdbID1}{$chainID1.$chainID2}{'evalue2'} = $line[10];
         $result{$pdbID1}{$chainID1.$chainID2}{'length2'} = $line[11];
         $result{$pdbID1}{$chainID1.$chainID2}{'bitScore2'} = $line[12];
         $result{$pdbID1}{$chainID1.$chainID2}{'organismMT2'} = $line[13];
         $c++;
      }
   close(FH);

   return(\%result);

}

#little Helper to get ids from dunbrack Fasta >1RYPA


sub parseForPdbAndChainID($){
   my $class = shift @_;
   my $INfile = shift(@_);

   my %result;
   open (INFILE, "$INfile") or die "Can not open/access result file: '$INfile'\n$!\n";

      while (my $line = <INFILE>){
         if ($line =~ /^>\w{5}\s/o){
            my $pdbID = substr($line,1,4);
            my $chainID = substr($line,5,1);
            $result{$pdbID}{$chainID}++;
         }
         elsif($line =~ /^\w{5}\s/o){
            my $pdbID = substr($line,0,4);
            my $chainID = substr($line,4,1);
            $result{$pdbID}{$chainID}++;
         }
else{next}
      }

   close(INFILE);
   return(\%result);
}



#1IGY  A  D  25.8
sub createAtomSequence($$){
   my $class = shift @_;
   my $baseDir = shift @_;
   my $HRinput = shift @_;
#$HRinput{$pdbID}{$chainID}
# my $base = '/home/s0571283/dunbrack/Latest/';

   my (%chainSequence);

   foreach my $pdbID (keys %{$HRinput}){
      $pdbID =~ tr/[A-Z]/[a-z]/;
      my $filename = $baseDir . $pdbID . '.sc';

      open(SC,'<',$filename) or die "Cannot open/access '$filename'\n$!\n";
#   print STDERR "Opened: $filename\n"; 
         my $residue;
         while(my $line = <SC>){

            if($line =~ /^SEQCRD/){
               my $chain      = substr($line,7,1);
               if(not exists $HRinput->{$pdbID}->{$chain}){next;}
               my $missing    = substr($line,15,3);
               if( $missing   =~ /---/){$residue = 'X';}
               else{ $residue = substr($line,9,1);}
               $pdbID =~ tr/[a-z]/[A-Z]/;
               $chainSequence{$pdbID.$chain} .= $residue;
            }
         }
            close(SC);
         }
   return(\%chainSequence);
}
# IN:
# >d1dlwa_ a.1.1.1 (A:) Protozoan/bacterial hemoglobin {Ciliate (Paramecium caudatum) [TaxId: 5885]}
# slfeqlggqaavqavtaqfyaniqadatvatffngidmpnqtnktaaflcaalggpnawt
# grnlkevhanmgvsnaqfttvighlrsaltgagvaaalveqtvavaetvrgdvvtv

# OUT:
# d1urvb_' => {
#   'orgName' => 'Homo sapiens',
#   'sequence' => 'elseaerkavqamwarlyansedvgvailvrffvnfpsakqyfsqfkhmedplemerspqlrkhasrvmgalntvvenlhdpdkvssvlalvgkahalkhkvepvyfkilsgvilevvaeefasdfppetqrawaklrgliyshvtaaykevgw',
#   'taxid' => '9606',
#   'scopID' => 'a.1.1.2',
#   'fasta' => '>d1urvb_ a.1.1.2 (B:) Cytoglobin {Human (Homo sapiens) [TaxId: 9606]}
#   elseaerkavqamwarlyansedvgvailvrffvnfpsakqyfsqfkhmedplemerspqlrkhasrvmgalntvvenlhdpdkvssvlalvgkahalkhkvepvyfkilsgvilevvaeefasdfppetqrawaklrgliyshvtaaykevgw'
#}

sub parseAstral100($){
   my $class = shift;
   my $Lastral100 = shift @_;
   my %result;

   open(FH,"$Lastral100") or die "Can not open/access '$Lastral100'\n$!\n";
   my $readingSeq = 0;
   my (@line,$line);
      LINE:
      while($line = <FH>){
         SWITCH:{
            if (!$readingSeq && $line =~ /^>/){
               $readingSeq = 1;
           
               @line = split(/\s/,$line);#print Dumper(\@line);
               $line[0] =~ s/^>//o;
               $result{$line[0]}{'scopID'} = $line[1];
               $line =~ /TaxId:\s(\d*)\]/o;
               $result{$line[0]}{'taxid'} = $1;
           
               if( $line =~ /\{.*\((.*)\)/){
                  $result{$line[0]}{'orgName'} = $1;
               }
               elsif($line =~ /\{([^[]*)\[TaxId:/){
                  my $orgname = $1;
                  $orgname =~ s/\s+$//;
                  $result{$line[0]}{'orgName'} = $orgname;
               }
               $result{$line[0]}{'fasta'} = $line;
               $result{$line[0]}{'sequence'} = '';
               last SWITCH;
            };
            if ($readingSeq && $line !~ /^>/){
               chomp($line);
               $result{$line[0]}{'sequence'} .= $line;
               $result{$line[0]}{'fasta'} .= $line;
               last SWITCH;
            };
            if ($readingSeq && $line =~ /^>/){
               $readingSeq = 0;
               @line = ();
               redo LINE;
            };

         }
      }
      close(FH);
return(\%result);
}
#           0         1             2         3                 4         5       6 
# # Fields: query id, subject id, % identity, alignment length, s. start, s. end, evalue
# # 18 hits found
# P27958<<Genome_polyprotein>>[Hepatitis_C_virus_genotype_1a_(isolate_H)]{11108}   pdb|1CU1|B  91.89 629   17 645   0.0
# P27958<<Genome_polyprotein>>[Hepatitis_C_virus_genotype_1a_(isolate_H)]{11108}   pdb|1CU1|A  91.89 629   17 645   0.0

sub parseBlast7($){
   my $class = shift @_;
   my $LblastOut = shift @_;
   my %result;
   open(FH,"$LblastOut") or die "Can not open/access '$LblastOut'\n$!\n";
      while(my $line = <FH>){
         if($line =~ /^#/){next;}
         chomp($line);
         my @line = split(/\t/,$line);
         if (scalar(@line) != 7){die "$LblastOut has wrong format\n";}
         $result{$line[0]}{$line[1]}{'query'}   =$line[0];
         $result{$line[0]}{$line[1]}{'subject'} =$line[1];
         $result{$line[0]}{$line[1]}{'%id'}     =$line[2];
         $result{$line[0]}{$line[1]}{'alignmentLength'}=$line[3];
         $result{$line[0]}{$line[1]}{'sstart'}  =$line[4];
         $result{$line[0]}{$line[1]}{'send'}    =$line[5];
         $result{$line[0]}{$line[1]}{'Eval'}    = $line[6];
         $result{$line[0]}{$line[1]}{'slength'} = ($line[5] - $line[6])+1 ;
         
      }
   close(FH);
   return(\%result);
}
sub parseBlast7C($){
   my $class = shift @_;
   my $LblastOut = shift @_;
   my %result;
   my $c = 1;
   open(FH,"$LblastOut") or die "Can not open/access '$LblastOut'\n$!\n";
      while(my $line = <FH>){
         if($line =~ /^#/){next;}
         chomp($line);
         my @line = split(/\t/,$line);
         if (scalar(@line) != 7){die "$LblastOut has wrong format\n";}
         $result{$c}{'query'}=$line[0];
         $result{$c}{'subject'}=$line[1];
         $result{$c}{'%id'}=$line[2];
         $result{$c}{'alignmentLength'}=$line[3];
         $result{$c}{'sstart'}=$line[4];
         $result{$c}{'send'}=$line[5];
         $result{$c}{'Eval'}=$line[6];
         $c++;
      }
   close(FH);
   return(\%result);
}
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
#sp|P10323|ACRO_HUMAN 1BMLB 34.52 252   143   7  42 292   20 250   1e-38  155

sub parseUniprotLikeVsBisc($$){
   my ($class,$LblastOut) =  @_;
   my %result;
   my $c = 0;


   open(FH,"$LblastOut") or die "Can not open/access '$LblastOut'\n$!\n";
      while(my $line = <FH>){
         next if $line =~ /^#/;
         chomp($line);
         my ($upid,$db);
         my @line = split(/\t/,$line);#print Dumper(\@line);sleep(1); 
         if (scalar(@line) != 12){die "$LblastOut has wrong format\n";}
         
         if($line[0] =~ /^\w+\|\w+\|/){
            my @ids = split(/\|/,$line[0]);

            $upid  = $ids[1];    
            $db    = $ids[2];
         }
         elsif(($line[0] =~ /^(\w{6,9})/) && $line[0] !~ /\|/){
            $upid = $1;
            $db = 'unknown';
         }
             

         if($line[1] !~ /^\w{5}$/){print STDERR "Not a proper pdbID: $line[1]\n$line\n";die}
         my $pdbID   = substr($line[1],0,4);
         my $chainID = substr($line[1],4,1);
         my $pdbAndChainID = $pdbID . $chainID;
         
         $line[10] = uc($line[10]);#1e10 -> 1E-10

         $result{$upid}{$c}{$pdbAndChainID}{'query'}        = $upid;
         $result{$upid}{$c}{$pdbAndChainID}{'db'}           = $db;
         $result{$upid}{$c}{$pdbAndChainID}{'subject'}      = $line[1];
         $result{$upid}{$c}{$pdbAndChainID}{'%id'}          = $line[2];
         $result{$upid}{$c}{$pdbAndChainID}{'aliLength'}    = $line[3];
         $result{$upid}{$c}{$pdbAndChainID}{'mismatches'}   = $line[4];
         $result{$upid}{$c}{$pdbAndChainID}{'gapOpenings'}  = $line[5];
         $result{$upid}{$c}{$pdbAndChainID}{'q.start'}      = $line[6];
         $result{$upid}{$c}{$pdbAndChainID}{'q.end'}        = $line[7];
         $result{$upid}{$c}{$pdbAndChainID}{'s.start'}      = $line[8];
         $result{$upid}{$c}{$pdbAndChainID}{'s.end'}        = $line[9];
         $result{$upid}{$c}{$pdbAndChainID}{'evalue'}       = $line[10];
         $result{$upid}{$c}{$pdbAndChainID}{'bitScore'}     = $line[11];
         $result{$upid}{$c}{$pdbAndChainID}{'slength'}      = ($line[9] - $line[8])+1 ;
 #           if(($upid eq 'Q6H6T1') &&($pdbAndChainID eq '2NRUB')){print "$line\n"; }

         $c++;
      }
   close(FH);
#print Dumper($result{'Q6H6T1'});die; 
   return(\%result);
}
sub parsePlasmodiumVsBisc($$){
   my ($class,$LblastOut) =  @_;
   my %result;
   my $c = 0;


   open(FH,"$LblastOut") or die "Can not open/access '$LblastOut'\n$!\n";
      while(my $line = <FH>){
         next if $line =~ /^#/;
         chomp($line);
         my ($upid,$db);
         my @line = split(/\t/,$line);#print Dumper(\@line);sleep(1); 
         if (scalar(@line) != 12){die "$LblastOut has wrong format\n";}
         
         my $id = $line[0]; 

         my $pdbID   = substr($line[1],0,4);
         my $chainID = substr($line[1],4,1);
         my $pdbAndChainID = $pdbID . $chainID;
         
         $line[10] = uc($line[10]);#1e10 -> 1E-10

         $result{$id}{$c}{$pdbAndChainID}{'query'}        = $id;
         $result{$id}{$c}{$pdbAndChainID}{'subject'}      = $line[1];
         $result{$id}{$c}{$pdbAndChainID}{'%id'}          = $line[2];
         $result{$id}{$c}{$pdbAndChainID}{'aliLength'}    = $line[3];
         $result{$id}{$c}{$pdbAndChainID}{'mismatches'}   = $line[4];
         $result{$id}{$c}{$pdbAndChainID}{'gapOpenings'}  = $line[5];
         $result{$id}{$c}{$pdbAndChainID}{'q.start'}      = $line[6];
         $result{$id}{$c}{$pdbAndChainID}{'q.end'}        = $line[7];
         $result{$id}{$c}{$pdbAndChainID}{'s.start'}      = $line[8];
         $result{$id}{$c}{$pdbAndChainID}{'s.end'}        = $line[9];
         $result{$id}{$c}{$pdbAndChainID}{'evalue'}       = $line[10];
         $result{$id}{$c}{$pdbAndChainID}{'bitScore'}     = $line[11];
         $result{$id}{$c}{$pdbAndChainID}{'slength'}      = ($line[9] - $line[8])+1 ;
 #           if(($upid eq 'Q6H6T1') &&($pdbAndChainID eq '2NRUB')){print "$line\n"; }

         $c++;
      }
   close(FH);
#print Dumper($result{'Q6H6T1'});die; 
   return(\%result);
}
# BLASTP 2.2.19+
# Query:  11BAA 124 XRAY 2.06 0.184 NA no RIBONUCLEASE, SEMINAL <GB S81747> [BOS TAURUS]
# Database: db/astral100SeqRes
# Fields: query id, subject id, % identity, q. start, q. end, s. start, s. end
# 267 hits found
#11BAA   d1tq9b_ 100.00  1       124     1       124
#11BAA   d1tq9a_ 100.00  1       124     1       124

sub parseBiscVsAstralBlastOutput($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^#/);
         my @line = split(/\t/,$line);
         unless ($line[2] == 100.00){next;}
         my $pdbID = substr($line[0],0,4);
         my $chainID = substr($line[0],4,1);
         my $sname = $line[1];

         push(@{$result{$pdbID.$chainID}},$sname);
      }
      close(FH);
   return(\%result);
}

#Query,  Subject, %id      alignLength,   mismatch,   gapOpenings,   q.start, q.end,   s.start, s.end,   e.val,   bit    
#1IGYA,  d1igya1, 100.00,  106,           0,          0,             1,       106,     1,       106,     1e-58,   220
#query and subject can sometimes have tab separated values->more than 12 fields! Check!!
sub parseBlastM9($){
   my $class = shift @_;
   my $LblastOut = shift @_;
   my %result;
   my $c = 1;
   open(FH,"$LblastOut") or die "Can not open/access '$LblastOut'\n$!\n";
      while(my $line = <FH>){
         next if $line =~ /^#/;
         chomp($line);
         my @line = split(/\t/,$line);
         if (scalar(@line) != 12){die "$LblastOut has wrong format.Fields in line (instead of 12): " . scalar(@line) ;print Dumper(\@line);die }
         $line[10] =~ s/e/E/o;
         $result{$c}{'query'}             = $line[0];
         $result{$c}{'subject'}           = $line[1];
         $result{$c}{'%id'}               = $line[2];
         $result{$c}{'alignmentLength'}   = $line[3];
         $result{$c}{'mismatches'}        = $line[4];
         $result{$c}{'gapOpenings'}       = $line[5];
         $result{$c}{'q.start'}           = $line[6];
         $result{$c}{'q.end'}             = $line[7];
         $result{$c}{'s.start'}           = $line[8];
         $result{$c}{'s.end'}             = $line[9];
         $result{$c}{'Eval'}              = $line[10];
         $result{$c}{'bitScore'}          = $line[11];
         $c++;
      }
   close(FH);
   return(\%result);
}
# IN:
#          10        20        30        40        50        60        70
# 01345678901234567890123456789012345678901234567890123456789012345678901234567890
# RES LYS B1549   141.30  70.4 132.67  81.2   8.63  23.0  94.60  81.2  46.70  55.4
sub parseNaccessRsa($){
   my ($class, $Lfile) = @_;
   my %result;
   open(FH,"$Lfile") or warn "Can not open/access '$Lfile'\n$!";
      while(my $line = <FH>){
         chomp($line);
         if($line =~ /^RES/){

            my $residueName = substr($line,4,3); $residueName =~ s/\s+//;
            my $chainID = substr($line,8,1); $chainID =~ s/\s+//;
            my $residueNumber = substr($line,9,4); $residueNumber =~ s/\s+//;
            my $asa = substr($line,14,8); $asa =~ s/\s+//;

            $result{'RES'}{$chainID}{$residueNumber}{$residueName} = $asa;
         }
         elsif($line =~ /^CHAIN/){
            my @line = split(/\s+/,$line);
            my $chainID = $line[2];
            my $asa = $line[3];
            $result{'CHAIN'}{$chainID} = $asa;
         }
         elsif($line =~ /^TOTAL/){
            my @line = split(/\s+/,$line);
            $result{'TOTAL'} = $line[1];
         }
      }
   close(FH);
   return(\%result); 
   
}

#00 serial        '41',
#01 atomName      'CG',
#02 altLocation   ' ' 
#03 residueName   'GLU',
#04 chainID       '1'
#05 resNo         '5',
#06 Xcord         '-3.053',
#07 Ycord         '194.754',
#08 Zcord         '103.857',
#09 accessiblity  '21.662',
#10 radius        '1.87',
#ATOM     41  CG  GLU 1   5      -3.053 194.754 103.857  21.662  1.87 
#ATOM     42  CD  GLU 1   5      -2.938 193.238 103.831  11.547  1.76 
#ATOM     43  OE1 GLU 1   5      -3.874 192.583 103.302  13.130  1.40 
#ATOM     44  OE2 GLU 1   5      -1.912 192.714 104.339  41.369  1.40 
#ATOM     45  C   GLU 1   5      -5.571 195.842 105.094   0.000  1.76 
#ATOM     46  O   GLU 1   5      -5.610 197.064 104.963   0.000  1.40 

sub parseNaccessAsa($){
   my $class = shift @_;
   my $Lfile =shift @_;
   my @result;

   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!\n";
      while(my $line = <FH>){
         chomp($line);
         if ($line =~ /^ATOM/){
            my @temp;
            my $serial      = substr($line,7,4);  $serial      =~ s/^\s+|\s+$//g; #0
            my $atomName    = substr($line,13,3); $atomName    =~ s/^\s+|\s+$//g; #1
            my $altLocation = substr($line,16,1); $altLocation =~ s/^\s+|\s+$//g; #2
            my $residueName = substr($line,17,3); $residueName =~ s/^\s+|\s+$//g; #3
            my $chain       = substr($line,21,1); $chain       =~ s/^\s+|\s+$//g; #4
            my $resNo       = substr($line,23,3); $resNo       =~ s/^\s+|\s+$//g; #5
            my $x           = substr($line,31,7); $x           =~ s/^\s+|\s+$//g; #6
            my $y           = substr($line,39,7); $y           =~ s/^\s+|\s+$//g; #7
            my $z           = substr($line,47,7); $z           =~ s/^\s+|\s+$//g; #8
            my $accessiblity= substr($line,55,7); $accessiblity=~ s/^\s+|\s+$//g; #9
            my $vdw         = substr($line,63,5); $vdw         =~ s/^\s+|\s+$//g; #10

            push(@temp,$serial);
            push(@temp,$atomName);
            push(@temp,$altLocation);
            push(@temp,$residueName);
            push(@temp,$chain);
            push(@temp,$resNo);
            push(@temp,$x);
            push(@temp,$y);
            push(@temp,$z);
            push(@temp,$accessiblity);
            push(@temp,$vdw);
            push(@temp,$line);
            push(@result,\@temp);

         }
      }
   close(FH);
#print Dumper(\@result);
return(\@result);

}
#  HETATM    1  C   ACE A   1      72.700   7.366   1.277  1.00  0.00           C  
#  HETATM    2  O   ACE A   1      72.796   6.159   1.052  1.00  0.00           O  
#  HETATM    3  CH3 ACE A   1      71.887   7.821   2.479  1.00  0.00           C  
#  HETATM    4  H1  ACE A   1      71.424   6.959   2.956  1.00  0.00           H  
#  HETATM    5  H2  ACE A   1      71.108   8.510   2.154  1.00  0.00           H  
#  HETATM    6  H3  ACE A   1      72.543   8.321   3.190  1.00  0.00           H  
#  ATOM      7  N   PRO A   2      73.287   8.299   0.502  1.00  0.00           N  
#  ATOM      8  CA  PRO A   2      74.031   7.998  -0.728  1.00  0.00           C  
#  ATOM      9  C   PRO A   2      75.461   7.456  -0.515  1.00  0.00           C  
############# take care ######################

=begin comment
sub parseATOMseq($){
   my ($class,$Lfile) =  @_;  
   my (%seq);
   my ($model,$endMdl) = (0,0);

   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!";
   my @file = <FH>;
   close(FH);
#Read in the ATOM and HETATM records
   foreach my $line (@file){
      chomp($line);
      if ($line =~ /^MODEL/){$model++}
      if ($line =~ /^ENDMDL/){$endMdl++;}
      if ( ($line =~ /^ATOM/) || ($line =~ /^HETATM/) ){
# store records according to chain
         my $residueName = substr($line,17,3); $residueName =~ s/^\s+|\s+$//g; #3
         my $chain       = substr($line,21,1); $chain       =~ s/^\s+|\s+$//g; #4
#Create ATOM sequence by concatenating the 3-letter residue code
         $seq{$endMdl}{'seq'}{$chain} .= $residueName;
#Keep the original record by concatenating the lines         
         $seq{$endMdl}{'records'}{$chain}  .= "$line\n";
      }    
   }

#remove HETATM only chains
HETATM:
   foreach my $chain (sort keys %{$seq{'0'}{'seq'}}){
#separate lines 
      my @lines = split(/\n/,$seq{'0'}{'records'}{$chain});
      my $flag = 0; 
#Iterate through lines anc check if at least one line of the chain contains an ATOM record
      for my $line(@lines){
         if ($line =~ /^ATOM/){$flag = 1;}
      }    
#only HETATM records
      if ($flag == 0){
         my $seq = $seq{'0'}{'seq'}{$chain};
#         print "Removed $chain\n"; 
#         print Dumper(\@lines); 
         my $del1 = delete ($seq{'0'}{'records'}{$chain});
         my $del2 = delete ($seq{'0'}{'seq'}{$chain});
#      my $del3 = delete ($seq{'0'}{'chains'}{$seq});
#         print "$del1\n"; 
#         print "$del2\n"; 
         next HETATM;
      }    
#      elsif($flag == 1){
#         $flag = 0;
#      }
   }
   foreach my $chain (sort keys %{$seq{'0'}{'seq'}}){
      $seq{'0'}{'list'}{$chain}++;
      my $seq = $seq{'0'}{'seq'}{$chain};
      push(@{$seq{'0'}{'chains'}{$seq}},$chain);
   }
#   print Dumper($seq{'0'}); 
   return($seq{'0'});
}

=cut





#  HETATM    1  C   ACE A   1      72.700   7.366   1.277  1.00  0.00           C  
#  HETATM    2  O   ACE A   1      72.796   6.159   1.052  1.00  0.00           O  
#  HETATM    3  CH3 ACE A   1      71.887   7.821   2.479  1.00  0.00           C  
#  HETATM    4  H1  ACE A   1      71.424   6.959   2.956  1.00  0.00           H  
#  HETATM    5  H2  ACE A   1      71.108   8.510   2.154  1.00  0.00           H  
#  HETATM    6  H3  ACE A   1      72.543   8.321   3.190  1.00  0.00           H  
#  ATOM      7  N   PRO A   2      73.287   8.299   0.502  1.00  0.00           N  
#  ATOM      8  CA  PRO A   2      74.031   7.998  -0.728  1.00  0.00           C  
#  ATOM      9  C   PRO A   2      75.461   7.456  -0.515  1.00  0.00           C  
############# take care ######################
sub parsePDBatom($){
   my ($class,$Lfile) =  @_;
#Translation of all 3-letter AA to 1-letter-codes used by Roland Dunbrack in his S2C service 
   my %letter =(  
         'NIY','Y','TPO','T','PCA','Q','PRO','P','CSS','C','YOF','Y','OCS','C','ILE','I','MEN','N','LLY','K','NEP','H','CGU','E','CYS','C','CSP','C','HIC','H','SCY','C','GLN','Q','CYG','C','ASN','N','VAL','V','KCX','K','TYR','Y','SME','M','GLX','Q','CMT','C','GL3','G','MLY','K','ASP','D','CSB','C','MHO','M','TYY','Y','TYI','Y','CSO','C','PHE','F','OCY','C','OPR','R','MGN','Q','CAS','C','LLP','K','MSE','M','PRS','P','PTH','Y','MHS','H','THR','T','AGM','R','TRN','W','SMC','C','CME','C','TYN','Y','ASX','N','GLY','G','HIP','H','PHD','D','CSE','C','TRF','W','IAS','D','SNC','C','SER','S','FME','M','GLU','E','HTR','W','FTR','W','ALA','A','MET','M','TYS','Y','MIS','S','NPH','C','HIS','H','MLZ','K','CSW','C','OMT','M','SVA','S','AYA','A','FGL','C','LYS','K','OLT','T','HYP','P','ALY','K','TRP','W','ACY','G','M3L','K','PYX','C','TYQ','Y','TPQ','Y','LEU','L','AAR','R','ARG','R','TRO','W','BHD','D','CSX','C','LYZ','K','CXM','M','AEI','T','PAQ','Y','CSD','C','SEP','S' 
         );
   my (%result,%altRes);
   my ($model,$endMdl) = (0,0);

   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!";
      my @file = <FH>;

   close(FH);

   foreach my $line (@file){
         chomp($line);

         if ($line =~ /^MODEL/){$model++}
         if ($line =~ /^ENDMDL/){$endMdl++;}
         if ($line =~ /^HETATM/){
            my $residueName = substr($line,17,3); $residueName =~ s/^\s+|\s+$//g; #3
               if($residueName ne 'MSE'){
                  my $chain       = substr($line,21,1); $chain       =~ s/^\s+|\s+$//g; #4
                  $altRes{$chain}{$residueName}++;# .= "$line\n";
                  next;
               }
         }
         if ( ($line =~ /^ATOM/) || ($line =~ /^HETATM/) ){

            my $serial      = substr($line,7,4);  $serial      =~ s/^\s+|\s+$//g; #0
            my $atomName    = substr($line,13,3); $atomName    =~ s/^\s+|\s+$//g; #1
            my $altLocation = substr($line,16,1); $altLocation =~ s/^\s+|\s+$//g; #2
            my $residueName = substr($line,17,3); $residueName =~ s/^\s+|\s+$//g; #3
            my $chain       = substr($line,21,1); $chain       =~ s/^\s+|\s+$//g; #4
            my $resNo       = substr($line,23,3); $resNo       =~ s/^\s+|\s+$//g; #5
            my $x           = substr($line,31,7); $x           =~ s/^\s+|\s+$//g; #6
            my $y           = substr($line,39,7); $y           =~ s/^\s+|\s+$//g; #7
            my $z           = substr($line,47,7); $z           =~ s/^\s+|\s+$//g; #8
            my $occupancy   = substr($line,55,5); $occupancy   =~ s/^\s+|\s+$//g; #9
            my $tempFactor  = substr($line,61,5); $tempFactor  =~ s/^\s+|\s+$//g; #10
            my $element     = substr($line,77,1); $element     =~ s/^\s+|\s+$//g; #11
#        my $charge      = substr($line,79,1); $charge      =~ s/^\s+|\s+$//g; #12
            $result{$endMdl}{$chain}  .= "$line\n";
         }
      }

   if($endMdl != $model){warn "<sarcasm>Surprisingly, something seems to be wrong with the PDB file $Lfile</sarcasm>\nInconsistensy between MODEL($model times) and ENDMODEL($endMdl times)";return}

   return(\%result,\%altRes);
}
sub parseMODELLERatom($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;
   my ($model,$endMdl) = (0,0);

   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!";
      my @file = <FH>;
   close(FH);

   foreach my $line (@file){
         chomp($line);
         if ($line =~ /^ATOM /){
            
            my $serial      = substr($line,7,4);  $serial      =~ s/^\s+|\s+$//g; #0
            my $atomName    = substr($line,13,3); $atomName    =~ s/^\s+|\s+$//g; #1
            my $altLocation = substr($line,16,1); $altLocation =~ s/^\s+|\s+$//g; #2
            my $residueName = substr($line,17,3); $residueName =~ s/^\s+|\s+$//g; #3
            my $chain       = substr($line,21,1); $chain       =~ s/^\s+|\s+$//g; #4
            my $resNo       = substr($line,23,3); $resNo       =~ s/^\s+|\s+$//g; #5
            my $x           = substr($line,31,7); $x           =~ s/^\s+|\s+$//g; #6
            my $y           = substr($line,39,7); $y           =~ s/^\s+|\s+$//g; #7
            my $z           = substr($line,47,7); $z           =~ s/^\s+|\s+$//g; #8
            my $occupancy   = substr($line,55,5); $occupancy   =~ s/^\s+|\s+$//g; #9
            my $tempFactor  = substr($line,61,5); $tempFactor  =~ s/^\s+|\s+$//g; #10
            my $element     = substr($line,77,1); $element     =~ s/^\s+|\s+$//g; #11
            $result{$chain}  .= "$line\n";
         }
      }
    return(\%result);
}
=begin comment
1  1  -  6        Record name   "ATOM  "
2  7  - 11        Integer       serial       Atom  serial number.
3  13 - 16        Atom          name         Atom name.
4  17             Character     altLoc       Alternate location indicator.
5  18 - 20        Residue name  resName      Residue name.
6  22             Character     chainID      Chain identifier.
7  23 - 26        Integer       resSeq       Residue sequence number.
8  27             AChar         iCode        Code for insertion of residues.
9  31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
10 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
11 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
12 55 - 60        Real(6.2)     occupancy    Occupancy.
13 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
14 77 - 78        LString(2)    element      Element symbol, right-justified.
15 79 - 80        LString(2)    charge       Charge  on the atom.
=cut
sub parsePDBatomResidue($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;
   my ($model,$endMdl) = (0,0);
   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!";
      my @file = <FH>;
   close(FH);

   my $Cline = 0;
   foreach my $line (@file){
      $Cline++;
      chomp($line);
      if ($line =~ /^MODEL/){$model++}
      if ($line =~ /^ENDMDL/){$endMdl++;}
      if($line =~  /^(ATOM  )(.....) (....)(.)(...) (.)(....)(.)   (........)(........)(........)(......)(......)           (..)/){
#print "1: $1 2: $2 3: $3 4: $4 5: $5 6: $6 7: $7 8: $8 9: $9 10: $10 11: $11 12: $12 13: $13 14: $14\n"; 

         my $iCode    = $8; 
         my $chain    = $6;
         my $resSeq   = $7;
         my $atomName = $3;

         $iCode    =~ s/^\s+|\s+$//g;
         $chain    =~ s/^\s+|\s+$//g;
         $resSeq   =~ s/^\s+|\s+$//g;
         $atomName =~ s/^\s+|\s+$//g;

            if($iCode !~ /\w/){$iCode = '*'}  
#         if($resSeq < 1){print "$Lfile\t$Cline\t$resSeq\n\n"; next;};
         if (not defined $resSeq){print $Lfile;die;} 
         
         $result{$endMdl}{$chain}{$resSeq}{$iCode}{$atomName}  .= "$line\n";
      }
   }
   if($endMdl != $model){warn "<sarcasm>Surprisingly, something seems to be wrong with the PDB file $Lfile</sarcasm>\nInconsistensy between MODEL($model times) and ENDMODEL($endMdl times)";return}
   return(\%result);
}

sub parsePDB($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my @result;

   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!";
      while(my $line = <FH>){
         chomp($line);
         if ($line =~ /^ATOM /){
            my @temp;
            my $serial      = substr($line,7,4);  $serial      =~ s/^\s+|\s+$//g; #0
            my $atomName    = substr($line,13,3); $atomName    =~ s/^\s+|\s+$//g; #1
            my $altLocation = substr($line,16,1); $altLocation =~ s/^\s+|\s+$//g; #2
            my $residueName = substr($line,17,3); $residueName =~ s/^\s+|\s+$//g; #3
            my $chain       = substr($line,21,1); $chain       =~ s/^\s+|\s+$//g; #4
            my $resNo       = substr($line,23,3); $resNo       =~ s/^\s+|\s+$//g; #5
            my $x           = substr($line,31,7); $x           =~ s/^\s+|\s+$//g; #6
            my $y           = substr($line,39,7); $y           =~ s/^\s+|\s+$//g; #7
            my $z           = substr($line,47,7); $z           =~ s/^\s+|\s+$//g; #8
            my $occupancy   = substr($line,55,5); $occupancy   =~ s/^\s+|\s+$//g; #9
            my $tempFactor  = substr($line,61,5); $tempFactor  =~ s/^\s+|\s+$//g; #10
            my $element     = substr($line,77,1); $element     =~ s/^\s+|\s+$//g; #11
            my $charge      = substr($line,79,1); $charge      =~ s/^\s+|\s+$//g; #12

            push(@temp,$serial);
            push(@temp,$atomName);
            push(@temp,$altLocation);
            push(@temp,$residueName);
            push(@temp,$chain);
            push(@temp,$resNo);
            push(@temp,$x);
            push(@temp,$y);
            push(@temp,$z);
            push(@temp,$occupancy);
            push(@temp,$tempFactor);
            push(@temp,$element);
            push(@temp,$charge);
            push(@temp,$line);

            push(@result,\@temp);

         }
      }
   close(FH);
#print Dumper(\@result);
return(\@result);

}
#00 name          'ATOM',
#01 serial        '4769',
#02 atomName      'HB3',
#03 residueName   'LEU',
#04 resNo         '301',
#05 Xcord         '1.814',
#06 Ycord         '210.952',
#07 Zcord         '112.676',
#08 charge        '0.0974',
#09 radius        '1.4870',
#10 whole line    'ATOM   4769  HB3 LEU   301       1.814 210.952 112.676  0.0974 1.4870'
sub parsePQR($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my @result;
   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!\n";
      while(my $line = <FH>){
         chomp($line);
         if ($line =~ /^ATOM /){
            my @line = split(/\s+/,$line);
            push(@line,$line);
            push(@result,\@line);
         }
      } 
   close(FH);
   return(\@result);
}
#d1dlwa_ 1dlw    A:      a.1.1.1 14982   cl=46456,cf=46457,sf=46458,fa=46459,dm=46460,sp=46461,px=14982
sub parseDirClaScoTxt($){
   my $class = shift @_;
   my $LscopClas = shift @_;
   my %result;
   my @result;

   open(FH,"$LscopClas") or die "Can not open/access '$LscopClas'\n$!\n";
      while(my $line = <FH>){
         if($line =~ /^#/){next;}
         my @line = split(/\t/,$line);

         my $astralID = $line[0];
         my $pdbID = $line[1];
         my $chainID = $line[2];
            $chainID =~ s/:.*//;
         my $scopID = $line[3];
         $result{'astral'}{$astralID} =  $scopID;
         push(@{$result{'pdb'}{$pdbID.$chainID}}, $scopID);

      }   
   close(FH);
   return(\%result)
}

# IN:
# 12E8H b.1.1.2  b.1.1.1

# OUT:
#'1WVEB' => {
#   'd.58.32.1' => 1,
#   'd.145.1.1' => 1
#  },
#3C5ZC d1lnua1  b.1.1.2  Class II MHC alpha chain, C-terminal domain 
#3C5ZC d1muja1  b.1.1.2  Class II MHC alpha chain, C-terminal domain 
#3C5ZC d1muja2  d.19.1.1 Class II MHC alpha chain, N-terminal domain 
sub parseBiscToScop($){
   my $class = shift @_;
   my $Lbisc2Scop = shift @_;
   my %result;

   open(FH,"$Lbisc2Scop") or die "Can not open/access '$Lbisc2Scop'\n$!\n";
      while(my $line = <FH>){
         chomp($line);
         my @line = split(/\t/,$line);
         $result{'scopText'}{$line[0]} = $line[3];
         $result{'scopID'}{$line[0]}   = $line[2];
         $result{'astral'}{$line[0]}   = $line[1];
      }
   close(FH);
   return(\%result);
}

# IN:
# # dir.des.scop.txt 
# # SCOP release 1.73 (November 2007)  [File format version 1.00]
# # http://scop.mrc-lmb.cam.ac.uk/scop/
# # Copyright (c) 1994-2007 the scop authors; see http://scop.mrc-lmb.cam.ac.uk/scop/lic/copy.html
# 46456 cl a  -  All alpha proteins
# 46457 cf a.1   -  Globin-like
# 46458 sf a.1.1 -  Globin-like
# 46459 fa a.1.1.1  -  Truncated hemoglobin
# 46460 dm a.1.1.1  -  Protozoan/bacterial hemoglobin
# 46461 sp a.1.1.1  -  Ciliate (Paramecium caudatum) [TaxId: 5885]
# 14982 px a.1.1.1  d1dlwa_  1dlw A:

# OUT:
# 'c.66.1.11' => 'Type II DNA methylase',
# 'a.70.1.1' => 'N-terminal domain of the delta subunit of the F1F0-ATP synthase',
# 'j.3.1.2' => 'Insect defense peptide thanatin',

sub parseDirDesScopTxtNew($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my (%result);

   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!\n";
      while(my $line = <FH>){
         if($line =~ /^#/){next;}
         chomp($line);
         my @line = split(/\t/,$line);#print Dumper(\@line);sleep(1); 
         my $sunid = $line[0];
         my $level = $line[1];
         my $sccs  = $line[2];
         my $sname = $line[3];
         my $descr = $line[4];
         if($level eq 'cl'|| $level eq 'cf' || $level eq 'sf' || $level eq 'fa'  ){
            $result{'scop'}{$sunid}{'level'} = $level;
            $result{'scop'}{$sunid}{'sccs'}  = $sccs;
            $result{'scop'}{$sunid}{'descr'} = $descr;
         }
         if($level eq 'px'){
            push(@{$result{'domain'}{$sccs}},$sname);
         }
         
      }   
   close(FH);
   return(\%result);
   
}
sub parseDirDesScopTxt($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my (%result);

   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!\n";
      while(my $line = <FH>){
         if($line =~ /^#/){next;}
         chomp($line);
         my @line = split(/\t/,$line);
         if ($line[1] eq 'cl'){
            $result{'cl'}{$line[2]} = $line[4];
         }   
         elsif ($line[1] eq 'cf'){
            $result{'cf'}{$line[2]} = $line[4];
         }   
         elsif ($line[1] eq 'sf'){
            $result{'sf'}{$line[2]} = $line[4];
         }   
         elsif ($line[1] eq 'fa'){
            $result{'fa'}{$line[2]} = $line[4];
         }   
         else{next;}
      }   
   close(FH);
   return(\%result);
}

# IN:
# $line[0] = chain; $line[1] = residueNumber; $line[2] = residue Name;  $line[3] = rel.interface
# H  3  GLN   1.05

sub parseInterfaceFile($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;

   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!\n";
      while(my $line = <FH>){
         chomp($line);
         my @line = split(/\t/,$line);
         $result{$line[0]}{$line[1]}{$line[2]} = $line[3];
         $result{'jmolString'}{$line[0]} .= $line[2] . $line[1] . ':' . $line[0] . ', ';
      }
   close(FH);
   foreach my $key(keys %{$result{'jmolString'}}){
      $result{'jmolString'}{$key} =~ s/, $//;
   }
   
   return(\%result);

}

sub parseInterface($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;

   open(FH,"$Lfile") or die "Can not open/access '$Lfile'\n$!\n";
      while(my $line = <FH>){
         chomp($line);
         my @line = split(/\t/,$line);
         my @id = split(/_/,$line[0]);
         my @chains = split(//,$id[1]);
         $result{$id[0]}{$chains[0]}{$chains[1]} = $line[1];
         $result{$id[0]}{$chains[1]}{$chains[0]} = $line[1];
      }
      close(FH);

      return(\%result);

}
sub parseUniprotText($){
   my ($class,$Lfile) =  @_;
   my %result;
   my ($reading_entry,$foundAC,$foundDE,$foundOS,$foundSQ) = (0,0,0,0,0);
   my $c = 0;
   my $d = 0;
   my ($firstAcc,@altAccs,$description,$org,$seq,$taxID);

   open( FH, "$Lfile" ) or die "Could not find the FASTA file $Lfile!\n";
   my ($line,$id);
LINE:
   while(  $line = <FH> ){
      chomp( $line );
SWITCH: {

           if ( !$reading_entry && $line =~ /^ID/ ) {
              $reading_entry = 1;
              last SWITCH;
           };

           if ( $reading_entry && !$foundAC && $line =~ /^AC/ ) {
              $foundAC = 1;
              $d++;
              $line =~ s/^AC//;
              $line =~ s/\s//g;
              my @line = split(/;/,$line);

              if(scalar (@line) > 1 ){$c++};
              $firstAcc = shift (@line);
              if(scalar (@line) > 0 ){
                 $result{'altAccs'}{$firstAcc} = \@line;
                 for my $altAcc (@line){
                    $result{'mapAltAccs'}{$altAcc} = $firstAcc;
                 }
              }
              else{$result{'altAccs'}{$firstAcc} = 0}
              last SWITCH;
           };
           if ( $reading_entry && !$foundDE && $line =~ /^DE/ ) {
               $foundDE = 1;
               $line =~ /Full=(.*)$/;
               $description = $1;
               $description =~ s/;$//;
              last SWITCH;
           }

           if ( $reading_entry && $foundOS == 0 && $line =~ /^OS/ ) {
              $foundOS = 1;
              $line =~ /^OS\s+(.*)[\.\(]$/;
              if(!$1){
                 $line =~ /^OS\s+(.*)$/;
                 $org = $1;
                 $line = <FH>;
                 $line =~ /^OS\s+(.*)[\.\(]$/;
                 $org .= " $1";
              }
              else{
                 $org = $1;
                 if($org =~ /\./){$org =~ s/\.//;}
                 if($org =~ /\(/){$org =~ s/\(.*\)//g;}
                 if($org =~ /^\s|\s+$/){$org =~ s/^\s+|\s+$//g;}
              }
              last SWITCH;
           }

           if ($reading_entry && $line =~ /^OX/){
              if($line =~ /NCBI_TaxID=(\d*);/){
                 $taxID = $1;
              }
              last SWITCH;
           }


           if ( $line =~ /^\/\// ) {

              $org =~ s/ /_/g;
           $result{'taxid'}{$firstAcc} = $taxID;
           $result{'description'}{$firstAcc} = $description;
           $result{'organism'}{$firstAcc} = $org;
           $result{'id'}{$firstAcc}++;

           ($reading_entry,$foundAC,$foundDE,$foundOS,$foundSQ) = (0,0,0,0,0);
           ($firstAcc,$description,$org,$seq,$taxID) = (0,0,0,0,0);
           undef(@altAccs);
           print "."; 

        };
        }
   }
   close(FH);
   print "$c/$d\n";
   return(\%result);
}


#AC   Q16653; O00713; O00714; O00715; Q13054; Q13055; Q14855; Q92891;
#AC   Q92892; Q92893; Q92894; Q92895; Q93053; Q96KU9; Q96KV0; Q96KV1;

sub parseFastaOrdered($){
   my $class = shift @_;
   my $LfastaFile = shift@_;
   my %result;
   my $reading_seq = 0;
   my $counter = 0;

   open( SOURCE, "$LfastaFile" ) or die "Could not find the FAASTA file $LfastaFile!\n";
   my ($line,$id);
LINE: 
   while(  $line = <SOURCE> ){
# print $line;sleep(1);
      chomp( $line );
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1;
              $id = $line;
              $id =~ s/^>//;

              $result{$counter}{'sequence'} = '';
              $result{$counter}{'fasta'} = "$line\n";
              $result{$counter}{'commentLine'} = $line;
              $result{$counter}{'id'} = $id;

              last SWITCH;
           };

           if ( $reading_seq && $line !~ /^>/ ) {

              $line =~ s/\s+//g;
              $result{$counter}{'sequence'} .= $line;
              $result{$counter}{'fasta'}    .= $line;
              last SWITCH;
           };

           if ( $reading_seq && $line =~ /^>/ ) {
              $reading_seq = 0;
              $counter++;
              redo LINE;
           };
        }
   }
   close(SOURCE);
   return (\%result);
}
sub parseFasta($){
   my $class = shift @_;
   my $LfastaFile = shift@_;
   my %result;
   my $reading_seq = 0;

   open( SOURCE, "$LfastaFile" ) or die "Could not find the FAASTA file $LfastaFile!\n";
   my ($line,$id);
LINE: 
   while(  $line = <SOURCE> ){
# print $line;sleep(1);
      chomp( $line );
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1;
              $id = $line;
              $id =~ s/^>//;

              $result{'sequence'}{$id} = '';
              $result{'fasta'}{$id} = "$line\n";
              $result{'commentLine'}{$id} = $line;

              last SWITCH;
           };

           if ( $reading_seq && $line !~ /^>/ ) {

              $line =~ s/\s+//g;
              $result{'sequence'}{$id} .= $line;
              $result{'fasta'}{$id} .= $line;
              last SWITCH;
           };

           if ( $reading_seq && $line =~ /^>/ ) {
              $reading_seq = 0;
              redo LINE;
           };
        }
   }
   close(SOURCE);
   return (\%result);
}
#>d1dlwa_ a.1.1.1 (A:) Protozoan/bacterial hemoglobin {Ciliate (Paramecium caudatum) [TaxId: 5885]}
#slfeqlggqaavqavtaqfyaniqadatvatffngidmpnqtnktaaflcaalggpnawt
#grnlkevhanmgvsnaqfttvighlrsaltgagvaaalveqtvavaetvrgdvvtv
#>d1dlya_ a.1.1.1 (A:) Protozoan/bacterial hemoglobin {Green alga (Chlamydomonas eugametos) [TaxId: 3053]}
#slfaklggreaveaavdkfynkivadptvstyfsntdmkvqrskqfaflayalggasewk
#gkdmrtahkdlvphlsdvhfqavarhlsdtltelgvppeditdamavvastrtevlnmpq
#q
#>d1s69a_ a.1.1.1 (A:) Protozoan/bacterial hemoglobin {Cyanobacteria (Synechocystis sp.), pcc 6803 [TaxId: 1143]}
#stlyeklggttavdlavdkfyervlqddrikhffadvdmakqrahqkafltyafggtdky
#dgrymreahkelvenhglngehfdavaedllatlkemgvpedliaevaavagapahkrdv
#
sub parseAstralFasta($){
   my ($class,$LfastaFile) = @_;
   my (%result,$line,$sname,$sccs);
   my $reading_seq = 0; 

   open(FH,'<', $LfastaFile ) or die "Could not find the FASTA file $LfastaFile!\n";
LINE: 
   while(  $line = <FH> ){
# print $line;sleep(1);
      chomp( $line );
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1; 
              my @line = split(/\s/,$line);
              $sname = $line[0];#d1dlwa_
              $sccs = $line[1];#a.1.1.1

              $result{'sequence'}{$sname} = '';
              $result{'sccs'}{$sname} = $sccs;
              $result{'fasta'}{$sname} = "$line\n";
              $result{'commentLine'}{$sname} = $line;

              last SWITCH;
           };   

           if ( $reading_seq && $line !~ /^>/ ) {

              $line =~ s/\s+//g;
              $result{'sequence'}{$sname} .= $line;
              $result{'fasta'}{$sname} .= $line;
              last SWITCH;
           };   

           if ( $reading_seq && $line =~ /^>/ ) {
              $reading_seq = 0; 
              $sname = '';
              $sccs = '';
              redo LINE;
           };   
        }    
   }
   close(FH);
   return (\%result);
}
sub parseUniprotFasta($){
   my ($class, $LfastaFile) = @_;
   my %result;
   my $reading_seq = 0; 

   open(FH,'<', $LfastaFile ) or die "Could not find the FASTA file $LfastaFile!\n";
   my ($line,$id);
LINE: 
   while(  $line = <FH> ){
# print $line;sleep(1);
      chomp( $line );
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1; 
              my @line = split(/\|/,$line);
              $id = $line[1];

              $result{'sequence'}{$id} = '';
              $result{'fasta'}{$id} = "$line\n";
              $result{'commentLine'}{$id} = $line;

              last SWITCH;
           };   

           if ( $reading_seq && $line !~ /^>/ ) {

              $line =~ s/\s+//g;
              $result{'sequence'}{$id} .= $line;
              $result{'fasta'}{$id} .= $line;
              last SWITCH;
           };   

           if ( $reading_seq && $line =~ /^>/ ) {
              $reading_seq = 0; 
              redo LINE;
           };   
        }    
   }
   close(FH);
   return (\%result);
}
sub parseUniprotFastaMultiple($$){
   my ($class, $LfastaFile,$HRresult) = @_;
   my $reading_seq = 0; 

   open(FH,'<', $LfastaFile ) or die "Could not find the FASTA file $LfastaFile!\n";
   my ($line,$id);
LINE: 
   while(  $line = <FH> ){
# print $line;sleep(1);
      chomp( $line );
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1; 
              $line =~ /^>.+\|(.+)\|/;
#my @line = split(/\|/,$line);
#$id = $line[1];
              $id = $1;

              $HRresult->{'sequence'}->{$id} = '';
              $HRresult->{'fasta'}->{$id} = "$line\n";
              $HRresult->{'commentLine'}->{$id} = $line;

              last SWITCH;
           };   

           if ( $reading_seq && $line !~ /^>/ ) {

              $line =~ s/\s+//g;
              $HRresult->{'sequence'}->{$id} .= $line;
              $HRresult->{'fasta'}->{$id} .= $line;
              last SWITCH;
           };   

           if ( $reading_seq && $line =~ /^>/ ) {
              $reading_seq = 0; 
              redo LINE;
           };   
        }    
   }
   close(FH);
}

sub parseHPRDfasta($){
   my $class = shift @_;
   my $LfastaFile = shift@_;
   my %result;
   my $reading_seq = 0;

   open( SOURCE, "$LfastaFile" ) or die "Could not find the FAASTA file $LfastaFile!\n";
   my ($line,$id);
LINE: 
   while(  $line = <SOURCE> ){
# print $line;sleep(1);
      chomp( $line );
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1;
              my @line = split(/\|/,$line);
              $id = $line[0];
              $id =~ s/^>//;

              $result{'sequence'}{$id} = '';
              $result{'fasta'}{$id} = "$line[0]\n";
              $result{'isoform'}{$id} = "$line[1]";

              last SWITCH;
           };

           if ( $reading_seq && $line !~ /^>/ ) {

              $line =~ s/\s+//g;
              $result{'sequence'}{$id} .= $line;
              $result{'fasta'}{$id} .= $line;
              last SWITCH;
           };

           if ( $reading_seq && $line =~ /^>/ ) {
              $reading_seq = 0;
              redo LINE;
           };
        }
   }
   close(SOURCE);
   return (\%result);
}
#>sp|Q65209|1002L_ASFB7 Protein MGF 100-2L OS=African swine fever virus (strain Badajoz 1971 Vero-adapted) GN=BA71V-151 PE=3 SV=1
#MGNKESKYLEMCSEEAWLNIPNIFKCIFIRKLFYNKWLKYQEKKLKKSLKLLSFYHPKKD

sub parseSwissProtFasta($){
   my $class = shift @_;
   my $LfastaFile = shift@_;
   my %result;
   my $reading_seq = 0;

   open( SOURCE, "$LfastaFile" ) or die "Could not find the FAASTA file $LfastaFile!\n";
   my ($line,$id);
LINE: 
   while(  $line = <SOURCE> ){
# print $line;sleep(1);
      chomp( $line );
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1;
#                        >sp|Q65209|1002L_ASFB7
              $line =~ /^>(\w*)\|(\w*)\|(\w*)_(\w*)\s/;

              $id = $2;
              my $org = $4;

              $result{'sequence'}{$id} = '';
              $result{'fasta'}{$id} = "$line\n";
              $result{'commentLine'}{$id} = $line;
              $result{'org'}{$id} = $org;

              last SWITCH;
           };

           if ( $reading_seq && $line !~ /^>/ ) {

              $line =~ s/\s+//g;
              $result{'sequence'}{$id} .= $line;
              $result{'fasta'}{$id} .= $line;
              last SWITCH;
           };

           if ( $reading_seq && $line =~ /^>/ ) {
              $reading_seq = 0;
              redo LINE;
           };
        }
   }
   close(SOURCE);
   return (\%result);
}
=cut
>uniprotkb:P25296/Calcineurin subunit B/Saccharomyces cerevisiae/4932
MGAAPSKIVDGLLEDTNFDRDEIERLRKRFMKLDRDSSGSIDKNEFMSIPGVSSNPLAGRIMEVFDADNSGDVDFQEFIT
GLSIFSGRGSKDEKLRFAFKIYDIDKDGFISNGELFIVLKIMVGSNLDDEQLQQIVDRTIVENDSDGDGRLSFEEFKNAI
ETTEVAKSLTLQYDV

>uniprotkb:O76061/Stanniocalcin-2/Homo sapiens/9606
MCAERLGQFMTLALVLATFDPARGTDATNPPEGPQDRSSQQKGRLSLQNTAEIQHCLVNAGDVGCGVFECFENNSCEIRG
LHGICMTFLHNAGKFDAQGKSFIKDALKCKAHALRHRFGCISRKCPAIREMVSQLQRECYLKHDLCAAAQENTRVIVEMI
HFKDLLLHEPYVDLVNLLLTCGEEVKEAITHSVQVQCEQNWGSLCSILSFCTSAIQKPPTAPPERQPQVDRTKLSRAHHG
EAGHHLPEPSSRETGRGAKGERGSKSHPNAHARGRVGGLGAQGPSGSSEWEDEQSEYSDIRR

=cut
sub parseIntActFasta($){
   my $class = shift @_;
   my $IntActFasta = shift @_;
   my %result;
   my $reading_seq = 0;

   open(FH,"$IntActFasta") or die "Can not open/access '$IntActFasta'\n$!\n";
   my ($line,$id);
LINE:
   while($line = <FH>){
      chomp($line);
SWITCH: {

           if ( !$reading_seq && $line =~ /^>/ ) {
              $reading_seq = 1;
              my @line = split(/\//,$line);
              $id = $line[0];
              $id =~ s/^>//;
              $result{'fasta'}{$id} = "$line\n";
              $result{'commentLine'}{$id} = $line;
              $result{'title'}{$id} = $line[1];
              $result{'organism'}{$id} = $line[2];
              $result{'taxid'}{$id} = $line[3];

              last SWITCH;
           };

           if ( ($reading_seq == 1) && ($line =~ /^[A-Z]/) ) {

              $line =~ s/\s+//g;
              $result{'sequence'}{$id} .= $line;
              $result{'fasta'}{$id} .= $line;
              last SWITCH;
           };

           if ( $reading_seq && $line =~ /^>/ ) {
              $reading_seq = 0;
              redo LINE;
           };
        }

   }
   close(FH);

   return(\%result);
}

# add the last sequence to the output list when
# we have run out of lines in the source file as
# this is the only one whose end is not delimited
# by the start of a new sequence
#
#      @entry = (
#            shift @header,
#            ( join ( " ", @header ) ),
#            ( join ( "", @sequence ) )
#            );
#      push @output, join ( "\t", @entry );
#
#
#      print ("Found " . @output . " FASTA files in $seq_source\n");
#
#      open ( DESTINATION, ">$seq_dest" ) or die
#         "Couldn't create the destination file $seq_dest!\n";
#
#      select ( DESTINATION );
#      foreach $entry ( @output ) {
#         print $entry . "\n";
#      }
#
#      close ( DESTINATION );

sub parseHB2file($){
   my $class = shift @_;
   my $Lhb2File = shift @_;
   my %result;
   my $c = 0;
   open(FH,"$Lhb2File") or die "Can not open/access '$Lhb2File'\n$!";
      while(my $line = <FH>){
         chomp($line);
         next if($line =~ /^HBPLUS/);
         next if($line =~ /^\(c\)/);
         next if($line =~ /^Citing/);
         next if($line =~ /^\s+</);
         next if($line =~ /^<-/);
         next if($line =~ /^c\s+i/);
         next if($line =~ /^h\s+n/);
         next if($line =~ /^n\s+s/);
$c++;
         my $chainA = substr($line,0,1); $chainA =~ s/^\s+|\s+$//g; #0
         my $resNoA = substr($line,1,4); $resNoA =~ s/^\s+|\s+$//g; $resNoA =~ s/0//g;#1
         my $insertionCodeA = substr($line,5,1); $insertionCodeA =~ s/^\s+|\s+$//g; #2
         my $resNameA = substr($line,6,3); $resNameA =~ s/^\s+|\s+$//g; #3
         my $atomNameA = substr($line,9,4);$atomNameA =~ s/^\s+|\s+$//g; #4

         my $chainB = substr($line,14,1); $chainB =~ s/^\s+|\s+$//g; #0
         my $resNoB = substr($line,15,4); $resNoB =~ s/^\s+|\s+$//g; $resNoB =~ s/0//g;#1
         my $insertionCodeB = substr($line,19,1); $insertionCodeB =~ s/^\s+|\s+$//g; #2
         my $resNameB = substr($line,20,3); $resNameB =~ s/^\s+|\s+$//g; #3
         my $atomNameB = substr($line,23,4);$atomNameB =~ s/^\s+|\s+$//g; #4

         $result{$c}{'chainA'} = $chainA;
         $result{$c}{'resNoA'} = $resNoA;
         $result{$c}{'insertionCodeA'} = $insertionCodeA;
         $result{$c}{'resNameA'} = $resNameA;
         $result{$c}{'atomNameA'} = $atomNameA; 

         $result{$c}{'chainB'} = $chainB;
         $result{$c}{'resNoB'} = $resNoB;
         $result{$c}{'insertionCodeB'} = $insertionCodeB;
         $result{$c}{'resNameB'} = $resNameB;
         $result{$c}{'atomNameB'} = $atomNameB; 

         $result{$c}{'line'} = $line; 


#print "$chainA\t$resNoA\t$insertionCodeA\t$resNameA\t$atomNameA\n";
#print "$chainB\t$resNoB\t$insertionCodeB\t$resNameB\t$atomNameB\n\n";


#print "$line\n";
      }
   close(FH);
   return(\%result);
}

#1  |  root  |     |  scientific name   |
#2  |  Bacteria |     |  scientific name   |
#6  |  Azorhizobium   |     |  scientific name   |
sub parseNcbiNamesDmp($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next unless($line =~ /scientific name/);
         chomp($line);
         my @line = split(/\|/,$line);
         my $taxid = $line[0];
         my $scientificName = $line[1];
         $taxid =~ s/^\s+|\s+$//g;
         $scientificName =~ s/^\s+|\s+$//g;
         $result{'taxid'}{$taxid} = $scientificName;
         $result{'scientificName'}{$scientificName} = $taxid;
      }   
   close(FH);
   return(\%result);
}
#code  |  name  |  preferred name |  taxid
#2  |  trimeresurus mucrosquamatus   |  Protobothrops mucrosquamatus  |  103944
#2  |  hordeum vulgare var. distichum   |  Hordeum vulgare subsp. vulgare   |  112509
# code 2 = taxid found, code 3 = not found;
#1 - the incoming name is our primary name for a taxon in our database
#2 - the incoming name is a secondary name for a taxon in our database
#(it could be listed as a synonym, a misspelling, a common name,
# or several other nametypes)
#3 - the incoming name is not found in our database
#+ - the incoming name is duplicated in our database
#(used in combination with the other status codes)

sub parseMissingBiscTaxids($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next unless($line =~ /^\d/);
         chomp($line);
         my @line = split(/\|/,$line);
         if(scalar (@line) != 4){die "Format of $Lfile wrong";}
         my $code          = $line[0];
         my $inputName     = $line[1];
         my $preferredName = $line[2];
         my $taxid         = $line[3];

         	$code          =~ s/^\s+|\s+$//g;
         	$inputName     =~ s/^\s+|\s+$//g;
         	$preferredName =~ s/^\s+|\s+$//g;
         	$taxid         =~ s/^\s+|\s+$//g;

            $inputName = lc($inputName);

         if(($code == 3) || ($code eq '+' )){warn "$code No taxid found for $inputName";next;}
         $result{$inputName}{'preferredName'} = $preferredName;
         $result{$inputName}{'taxid'} = $taxid;
      }
      return(\%result);
}

#uniprotkb:P25296  Saccharomyces cerevisiae   4932
#uniprotkb:P35381  Drosophila melanogaster 7227
#uniprotkb:Q0P7W7  Campylobacter jejuni 197
sub parseIntActToOrg($){
   my $class = shift @_;
   my $Lfile = shift @_; 
   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '" . $Lfile ."'\n$!";
      while(my $line = <FH>){
         chomp($line);
         my @line = split(/\t/,$line);
         $result{$line[0]}{'organism'} = $line[1];
         $result{$line[0]}{'taxid'} = $line[2];
      }   
   close(FH);
   return(\%result);
}
sub getPdbHeader($){
   my $class = shift @_;
   my $LpdbFile = shift @_; 
   my $header;
   open(FH,"$LpdbFile") or die "Can not open/access '$LpdbFile'\n$!\n";
   while(my $line = <FH>){
      if ($line =~ /^HEADER/){
         $header = substr($line,10,40);
         $header =~ s/^\s+|\s+$//g;
         last;
      }   
   }   
   close(FH);
   return($header);
}

sub getPdbTitle($){
   my $class = shift @_;
   my $LpdbFile = shift @_; 
   my $title;
   my $flag = 0;
   open(FH,"$LpdbFile") or die "Can not open/access '$LpdbFile'\n$!\n";
   while(my $line = <FH>){
      if (($line =~ /^TITLE/) && ($flag == 0)){
         my $string = substr($line,10,69);
         $string =~ s/\s+$//;
         $string =~ s/^\s+//;
         $title = $string ;
         $flag = 1;
      }   
      elsif(($line =~ /^TITLE/) && ($flag == 1) ){
         my $string = substr($line,11,39);
         $string =~ s/\s+$//;
         $string =~ s/^\s+//;
         $title .= " $string";
      }   
      elsif(($line !~ /^TITLE/) && ($flag == 1)){
         $flag = 0;
         last;
      }   
   }   
   close(FH);
   return($title);
}

sub executeNeedleWater($$$){
   my $class = shift @_;
   my ($SRseq1, $SRseq2, $SRembossDir, $SRalgorithm) = @_; 
   my $result;
   open(F1,">seq1") or die "Can not open/access 'seq1'\n$!\n" ;
   open(F2,">seq2") or die "Can not open/access 'seq2'\n$!\n" ;
   print F1 $$SRseq1;
   print F2 $$SRseq2;
   close(F1);
   close(F2);

   my $cmd = $$SRembossDir . $$SRalgorithm . " ";
   $cmd .= "seq1 seq2 -gapopen 10.0 -gapextend 0.5 -sprotein1 -sprotein2 -outfile temp.out";
   if ( system ($cmd) == 0){ print  "$$SRalgorithm had a happy end!\n"; }
   else{
      my $exitValue = $? / 256;
      print STDERR "Miserably failed to do something:$$SRalgorithm $exitValue\n\n";
      print "$cmd\n";
      die;
   }   

   open(FH,"<temp.out") or die "Can not open/access 'temp.out'\n$!\n";
LOOP:
   while(my $line = <FH>){
# Identity:    99/453 (21.9%)
      if($line =~ /^# Identity:/){
         $line =~ /\((.*)%\)$/;
         $result = $1; 
         last LOOP;
      }   
   }   
   close(FH);
   if ( system ("rm -f seq1 seq2 temp.out") == 0){ 
      print "removed seq1 seq2 temp.out\n";
   }   
   else {
      my $exitValue = $? / 256;
      print "Error removing files: $exitValue\n";
   }   
   return(\$result);

}

sub parseDir($$$){
   my ($class, $path, $pattern) =  @_;
   my @files;

   opendir(DIR,$path) or die 'Can not open dir ' . $path . "\n$!";
      @files = grep { /$pattern/   } readdir(DIR); 
   closedir(DIR);
   return(\@files);
}
sub stamp(){
   my $class = shift;
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,
         $yday,$isdst)= localtime(time);
   my $stamp =  sprintf("%4d-%02d-%02d|%02d:%02d:%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
   return($stamp);
}
sub stampDay(){
   my $class = shift;
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,
         $yday,$isdst)= localtime(time);
   my $stamp =  sprintf("%4d-%02d-%02d", $year+1900,$mon+1,$mday);
   return($stamp);
}


# 0       | 1        | 2                | 3            | 4               | 5                   | 6                    | 7                   | 8
# #PDB ID | \%seq id | alignment length | indentities  |  ASA interface  | ASA single chain 1  |  ASA complex chain 1 | ASA single chain 2  |  ASA complex chain 2
# 1a14:H:L 29 97 29 766   6216.2   5462.9   5515.0   4734.7
sub parseInteractions($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next unless($line =~ /^(\w{4}):(\w):(\w)/);
         chomp($line);
         my @line = split(/\t/,$line);
         my $pdbID = $1;
         my $chain1 = $2;
         my $chain2 = $3;
         my $id = $pdbID . ':'.$chain1 .':'. $chain2;
         $result{line}{$pdbID}{$chain1}{$chain2} = $line;
         $result{seqId}{$id}       = $line[1];
         $result{aliLength}{$id}   = $line[2];
         $result{indentities}{$id} = $line[3];
         $result{interface}{$id}   = $line[4];
         $result{single1}{$id}     = $line[5];
         $result{complex1}{$id}    = $line[6];
         $result{single2}{$id}     = $line[7];
         $result{complex2}{$id}    = $line[8];
         $result{ids}{$pdbID}{$chain1}++;
         $result{ids}{$pdbID}{$chain2}++;
         $result{file}{$id} = $pdbID.'_'.$chain1.$chain2;
         $result{dir}{$id} = lc(substr($pdbID,1,2)) . '/';
#         $result{dir}{$pdbID} = lc(substr($pdbID,1,2)) . '/';
         $result{pdbID}{$id} = $pdbID;
      }
   close(FH);

   return(\%result);
}
#1RJL  C  Borrelia burgdorferi 139
#2GTL  C  Lumbricus terrestris 639
sub parseBiscToOrg($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         if($line =~ /^#/){next;}
         chomp($line);
         my @line = split(/\t/,$line);
         my $pdbID = $line[0];
         my $chainID = $line[1];
         my $name = $line[2];
         my $taxid = $line[3];

         $result{$pdbID . $chainID}{'name'} = $name;
         $result{$pdbID . $chainID}{'taxid'} = $taxid;
      }
   close(FH);
   return(\%result);
}
# 3DGV  CRYSTAL STRUCTURE OF THROMBIN ACTIVATABLE FIBRINOLYSIS 2 INHIBITOR (TAFI)  HYDROLASE
# 1F4O  CRYSTAL STRUCTURE OF GRANCALCIN WITH BOUND CALCIUM METAL TRANSPORT

sub parseBiscToTitleHeader($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         if($line =~ /^#/){next;}
         chomp($line);
         my @line = split(/\t/,$line);
         my $pdbID  = $line[0];
         my $title  = $line[1];
         my $header = $line[2];
         $result{$pdbID}{'title'} = $title;
         $result{$pdbID}{'header'} = $header;
      }
   close(FH);
   return(\%result);

}

sub getAtomFasta{
   my $class = shift @_;
   my $Lfile = shift @_;

   my %result;
   my $lineno = 0;
   my %seq;

   my %aa_abbrev2letter = (
         "PRO" => "P",
         "GLY" => "G",
         "LEU" => "L",
         "LYS" => "K",
         "SER" => "S",
         "ALA" => "A",
         "PHE" => "F",
         "TRP" => "W",
         "GLN" => "Q",
         "GLU" => "E",
         "VAL" => "V",
         "ILE" => "I",
         "CYS" => "C",
         "TYR" => "Y",
         "HIS" => "H",
         "ARG" => "R",
         "ASN" => "N",
         "ASP" => "D",
         "THR" => "T",
         "MET" => "M",

         "A" => "A",
         "T" => "T",
         "C" => "C",
         "G" => "G",
         "U" => "U",
         );

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
   my $c = 0;
   while(my $line = <FH>){
      ++$lineno;
# get the sequence from the ATOM record
      if($line =~  /^(ATOM  )(.....) (....) (...) (.)(....)(.)/){
         my $altLoc  = $4;
         my $chainID = $5;
         my $resSeq  = $6;
         $seq{$chainID}[$0] = $aa_abbrev2letter{$altLoc} || "X";

#         if(! defined $seq{$chainID}->[$resSeq - 1]){
#            $seq{$chainID}->[$resSeq - 1]
#               = $aa_abbrev2letter{$altLoc} || "X";
#         }
#         elsif($seq{$chainID}->[$resSeq - 1] ne ($aa_abbrev2letter{$altLoc} || "X")){
#            die "$Lfile ATOM record at line $lineno indicates that residue " .
#               "at insertion position " . ($resSeq + 0) . " is $altLoc, " .
#               "which contradicts previous record.\n";
#         }
      }
   }
   foreach my $chainID (keys %seq){
      my $seqX = join "", map { $_ || "X" } @{$seq{$chainID}};
      $seqX =~ s/(.{80}|.{1,79}$)/$1\n/g;

      my $seqDash = join "", map { $_ || "-" } @{$seq{$chainID}};
      $seqDash =~ s/(.{80}|.{1,79}$)/$1\n/g;
      
      my $seqPlain = join "", map { $_ || '' } @{$seq{$chainID}};
      $seqPlain =~ s/(.{80}|.{1,79}$)/$1\n/g;

      $result{$chainID}{'X'} = $seqX;
      $result{$chainID}{'-'} = $seqDash;
      $result{$chainID}{'plain'} = $seqPlain;
   }

   return(\%result);
   close(FH);
}
# 1RYP<tab>A
sub parsePqrErrorFiles ($){
   my $class = shift @_;
   my $Lfile = shift @_;
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         my @line = split(/\t/,$line);
         $result{$line[0]}{$line[1]}++;
      }
   close(FH);
   
   return(\%result);

}
#HPRD_ID_MAPPINGS.txt
#hprd_id,geneSymbol,nucleotide_accession,protein_accession,entrezgene_id,omim_id,swissprot_id,main_name
#00001   ALDH1A1    NM_000689.3          NP_000680.2       216           0       P00352   Aldehyde dehydrogenase 1

sub parseHprdIDMap($){
   if(scalar @_ != 2){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 2 needed\n";die;}
   my ($class,$Lfile) = @_;

   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
   while(my $line = <FH>){
      chomp($line);
      my @line = split(/\t/,$line);
      my $id = $line[0];
      $result{$id}{'geneSymbol'} = $line[1];
      $result{$id}{'nucleotideAccession'} = $line[2];
      $result{$id}{'proteinAccession'} = $line[3];
      $result{$id}{'entrezgeneID'} = $line[4];
      $result{$id}{'omimID'} = $line[5];
      $result{$id}{'accession'} = $line[6]; #all other programs define accession as swissProt
      $result{$id}{'name'} = $line[7];

   }
   close(FH);
   return(\%result);
}
#BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt
#0                       1                    2                       3                       4                    5                      6               7
#interactor_1_geneSymbol,interactor_1_hprd_id,interactor_1_refseq_id, interactor_2_geneSymbol,interactor_2_hprd_id,interactor_2_refseq_id,experiment_type,         reference_id
#ALDH1A1,                00001	               NP_000680.2	            ALDH1A1	               00001	               NP_000680.2	           in vivo;yeast 2-hybrid	12081471,16189514
sub parseHprdPPI($){
   if(scalar @_ != 2){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 2 needed\n";die;}
   my ($class,$Lfile) = @_;
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
   while(my $line = <FH>){
      my @line = split($line); 
      $result{$line[1]}{$line[4]}++;
      $result{$line[4]}{$line[1]}++;
   }
   close(FH);
   return(\%result);
}
#  BIOGRID_ID      IDENTIFIER_VALUE        IDENTIFIER_TYPE ORGANISM_OFFICIAL_NAME
#  1       814566  ENTREZ_GENE     Arabidopsis thaliana
#  1       rrn26   OFFICIAL_SYMBOL Arabidopsis thaliana
#  1       ArthMr001       LOCUS_TAG       Arabidopsis thaliana
#  2       814568  ENTREZ_GENE     Arabidopsis thaliana
#  2       nad9    OFFICIAL_SYMBOL Arabidopsis thaliana
#  2       ArthMp007       LOCUS_TAG       Arabidopsis thaliana
#  2       CAA69753        PROTEIN_ACCESSION       Arabidopsis thaliana
#  2       6851007 PROTEIN_GI      Arabidopsis thaliana
#  2       NP_085479       PROTEIN_ACCESSION       Arabidopsis thaliana
#  2       13449296        PROTEIN_GI      Arabidopsis thaliana
#  2       1785679 PROTEIN_GI      Arabidopsis thaliana
#  2       Q95748  UNIPROTKB       Arabidopsis thaliana

#Multiple uniprot ids for one BioGrid possible
sub parseBioGridIDs($){
   if(scalar @_ != 2){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 1 needed\n";die;}
   my ($class,$Lfile) = @_;
   my (%result);
   my $flag = 0; 

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         if($line =~ /^(\d+)\t(\S+)\t(\S+)\t/){
            if(exists $result{'ids'}{$2}){die};
            $result{'biogrid'}{$2} = $1;
            $result{'type'}{$2} = $3;
            push(@{$result{'uniprot'}{$1}},$2) if ($3 eq 'UNIPROTKB');
         }
      }    
   close(FH);
   return(\%result);
}

#	INTERACTOR_A    INTERACTOR_B    OFFICIAL_SYMBOL_A       OFFICIAL_SYMBOL_B       ALIASES_FOR_A   ALIASES_FOR_B   EXPERIMENTAL_SYSTEM     SOURCE  PUBMED_ID
#	       ORGANISM_A_ID   ORGANISM_B_ID
#	EG6416  EG2318  MAP2K4  FLNC    JNKK|JNKK1|MAPKK4|MEK4|MKK4|PRKMK4|SEK1|SERK1   ABP-280|ABP280A|ABPA|ABPL|FLN2|FLJ10186 Two-hybrid      Marti A (1997)  90068
#	95 9606    9606
#	EG84665 EG88    MYPN    ACTN2   MYOP    N/A     Two-hybrid      Bang ML (2001)  11309420        9606    9606
sub parseBioGridPPI($){
   if(scalar @_ != 2){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 1 needed\n";die;}
   my ($class,$Lfile) = @_;
   my (%result);
   my $flag = 0; 
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         if($line =~ /^INTERACTOR_A/){$flag = 1;}
         if($flag == 1 and $line =~ /^(\S+)\t(\S+)/){
            $result{$1}{$2}++;
            $result{$2}{$1}++;
         }
      }
   close(FH);
   return(\%result);
}

sub parseParsedUniprot($){
   my ($class,$Lfile) = @_; 
   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
   while(my $line = <FH>){
      chomp($line);
      my @line = split(/\t/,$line); #[0]=primary acc,[1]=empty or alt acc, [2]=taxid
      $result{'acc'}{$line[0]}{'taxid'} = $line[2];
      if ( ref($line[1]) eq "ARRAY"){
         my @altAccs = split(/,/,$line[1]);
         for my $altAcc(@altAccs){
            $result{'altAcc'}{$altAcc}{$line[2]} = $line[0];
            $result{'mapAltAccs'}{$altAcc} = $line[0];
         }
      }   
   }   
   close(FH);
   return(\%result);
}

sub parseIntActinteractions($){
   my ($class,$Lfile) =  @_; 
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '" . $Lfile ."'\n$!";
      while(my $line = <FH>){
         chomp($line);
         my @line = split(/\t/,$line);
         for (my $index = 1 ; $index < scalar(@line) ; $index++){
            $result{$line[0]}{$line[$index]}++;
         }
      }   
   close(FH);
   return(\%result);
}

#get isoform ID (line[1]
sub parseHprdSeqInfo($){
   my ($class,$Lfile) =  @_;  
   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '" . $Lfile ."'\n$!";
      while(my $line = <FH>){
         chomp($line);
         my @line = split(/\t/,$line);
            $result{$line[0]}{$line[1]}++;
      }    
   close(FH);
   return(\%result);
}

sub parsePlasmosiumInteractions($){
   my ($class,$Lfile) =  @_; 
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
sub parsePlasmodiumInteractionsTitle($){
   my ($class,$Lfile) =  @_; 
   my %result;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^"bait/);
         chomp($line);
         $line =~ s/^"//;
         $line =~ s/"$//;
         my @line = split(/","/,$line);
         $result{$line[0]} = $line[1];
         $result{$line[2]} = $line[3];
      }   
   close(FH);
   return(\%result);
}
# PDB   CHAIN SP_PRIMARY  GO_ID
# 101m  A  P02185      GO:0015671
# Take care, 2 tabs between UniProtACC and GO
sub parsePdb2GO($){
   if((scalar (@_)-1 ) != 1){print STDERR "Argument list length: " . scalar(@_)-1 . " argument(s) passed, 1 needed\n";die;}
   my ($class,$Lfile) =  @_; 
   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^PDB/);
         chomp($line);
         my @line = split(/\t/,$line);
         $line[0] = uc($line[0]);
         $result{$line[0]}{$line[1]}{'uniprot'} = $line[2];
         push(@{$result{$line[0]}{$line[1]}{'GO'}},$line[4]);
      }
         
   return(\%result);
}

#11AS  A
sub parseBiscUniqueIDs($){
   if((scalar (@_)-1 ) != 1){print STDERR "Argument list length: " . scalar(@_)-1 . " argument(s) passed, 1 needed\n";die;}
   my ($class,$Lfile) =  @_; 
   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^PDB/);
         chomp($line);
         my @line = split(/\t/,$line);
         $result{$line[0]}{$line[1]}++;
      }
         
   return(\%result);
}
=begin comment
PISA v1.01 << 29/10/2007 >>   Session 1fnt

 Analysis of protein interfaces suggests that the following
 quaternary structures are stable in solution
 ----.-----.---------------------------------------.---------------
 Set |  No | Size  Id      ASA       BSA    DGdiss | Formula
 ----+-----+---------------------------------------+---------------
   1 |   1 |   7    0   69127.0   21818.3    106.0 | A(7)
     |   2 |   7    0   69142.0   21883.8    105.9 | A(7)
     |   3 |   7    1   71780.0   16984.8     11.9 | ABCDEFGa
     |   4 |   7    1   71809.0   16953.6     11.6 | ABCDEFGa
     |   5 |   4    2   34653.3   11520.2     16.2 | ABCDa(5)
     |   6 |   4    2   34712.5   11504.4     14.5 | ABCDa(5)
     |   7 |   2    3   19690.9    2708.1      2.0 | ABa
     |   8 |   2    3   19718.7    2679.0      1.5 | ABa
     |   9 |   1    4   10348.1       0.0     -0.0 | A
     |  10 |   1    4   10354.9       0.0     -0.0 | A
 ----+-----+---------------------------------------+---------------
=cut
sub parsePisaAssembly($){
   if((scalar (@_)-1 ) != 1){print STDERR "Argument list length: " . scalar(@_)-1 . " argument(s) passed, 1 needed\n";die;}
   my ($class,$Lfile) =  @_; 
   my (%result,$rec,$set,$stable);
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         if($line =~ /quaternary structures are stable in solution/){$stable = '1'}
         elsif($line =~ /quaternary structure is stable in solution/){$stable = '1'}
         elsif($line =~ /These structures may or may not be/){$stable = '0'}
         elsif($line =~ /This structure may or may not be/){$stable = '0'}
         elsif($line =~ /however they do not/){$stable = '-1'}
         elsif($line =~ /no calculation results found in session/){$stable = '-2';$result{$stable}{0}{0} = 'no calculation results found in session';next;}
         elsif($line =~ /complexate in solution/){$stable = '-3';$result{$stable}{0}{0} = 'Does not complexate in solution';next;}
         if($line =~ /^\s*(\d+)*\s*\|\s*(\d+)\s*\|\s*(\d+)\s*(\d)\s*(\S+)\s+(\S+)\s+(\S+)\s*\|\s*(\S+)\s+/){
            if($1){$set = $1}
            $rec = {
               Set => $set,
               Serial => $2,
               Size => $3,
               Id => $4,
               ASA => $5,
               BSA => $6,
               DGdiss => $7,
               Formula => $8,
               Stable => $stable
            };
            warn $Lfile if (not defined $rec->{Stable});
            warn $Lfile if (not defined $rec->{Set});
            warn $Lfile if (not defined $rec->{Size});
            next if(exists  $result{ $rec->{Stable} }{ $rec->{Set} }{ $rec->{Size} });
            $result{ $rec->{Stable} }{ $rec->{Set} }{ $rec->{Size} } = $rec;
         } 
      }
#      foreach my $set (sort {$a <=> $b} keys %result){
#         foreach my $serial (sort {$a <=> $b} keys %{$result{$set}}){
#            print "$set\t$serial\n"; 
#            print Dumper($result{$set}{$serial}); 
#         }
#      }
         
   return(\%result);

}
sub parseATOMseq(){
   if((scalar (@_)-1 ) != 1){print STDERR "Argument list length: " . scalar(@_)-1 . " argument(s) passed, 1 needed\n";die;}
   my ($class,$Lfile) =  @_; 
#   print "$Lfile\n"; 
#   $Lfile = '/home/tmp/BISC/update/pdb/101m.pdb';
#   $Lfile = '/home/tmp/BISC/update/pisa/structures/original/10gs.pisa.1.pdb';

   my (%finished,%temp,%test,@lines);
   my ($flag,$oldChain,$FfirstLine) = (0,-1,1); 

# Read assembly in array. Also X-ray files might have multiple models, only read in 1. Models are (or should be!) separated by ENDMDL record
   open( FH, '<' , $Lfile ) or die "Could not find $Lfile!\n";
      while(my $line = <FH> ){
         last if ($line =~ /^ENDMDL/ || $line =~ /^END/);
         next unless( ($line =~ /^ATOM/) || ($line =~ /^HETATM/) ||  ($line =~ /^TER/) );
         chomp( $line );
         push(@lines,$line);
      }
   close(FH);
#print "$lines[$#lines]\n"; 
LINE:
   for my $line( @lines){

#We are at the end of a chain
# Test if all entries are HETATM
      if($line =~ /^TER/ ){
            my ($Fredo) = &newChain($flag,\%finished,\$oldChain,\%temp,$Lfile,$line);
            if($Fredo == 1){redo LINE}
         $flag = 0;
         %temp = ();
         if($line eq $lines[$#lines]){
# We are at the end of the file, create a few more useful hashes
# Helpful for homooloigomers: Which chain IDs are assigned to a sequence              
            &end(%finished);
         }
         else{
            next LINE;
         }
      }#$line=TER
#Not TER      
      my $resName = substr($line,17,3); 
         $resName =~ s/^\s+|\s+$//g; #3
      my $newChain   = substr($line,21,1); 
      my $x          = substr($line,31,7); $x =~ s/^\s+|\s+$//g; #6
      my $y          = substr($line,39,7); $y =~ s/^\s+|\s+$//g; #7
      my $z          = substr($line,47,7); $z =~ s/^\s+|\s+$//g; #8
#Get the XYZ coordinates of the first atom of a new chain



#We are at the end.  In PISA files the last line still contains a record. The last line can either be a continuing chain or a chain on it s own (mostly a HETATM)
# The string comparison should work because the last line should contain X,Y & Z coordinates.
      if($line eq $lines[$#lines]) {
#         print "$Lfile\n"; 
#print "EOF\n"; 
# The last line can either continue a chain or be a chain on its own

# Continuation of chain
         if( $newChain eq $oldChain ){
            if($line =~ /^ATOM/){$flag = 1}

            $temp{$oldChain}{sequence} .= $resName;
            $temp{$oldChain}{records}  .= "$line\n";
            $temp{$oldChain}{coord}    .= $x . $y . $z;
            my ($Fredo) = &newChain($flag,\%finished,\$oldChain,\%temp,$Lfile,$line);
            if($Fredo == 1){redo LINE}

         }
# A new chain
         elsif( $newChain ne $oldChain  ){
            my ($Fredo) = &newChain($flag,\%finished,\$oldChain,\%temp,$Lfile,$line);
            if($Fredo == 1){redo LINE}
# First finish the old chain
# Start new record
            $oldChain = $newChain;
# Is the new chain a ATOM or HETATM record?
            if($line =~ /^ATOM/){$flag = 1}
            if($flag == 0){
               $finished{onlyHet}{$newChain}{sequence} = $resName;
               $finished{onlyHet}{$newChain}{records}  = $line;
               $temp{$newChain}{coord}                 = $x . $y . $z;
            }
            elsif($flag == 1){
               $finished{normal}{$newChain}{sequence} = $resName;
               $finished{normal}{$newChain}{records}  = $line;
               $temp{$newChain}{coord}                = $x . $y . $z;
               $finished{chains}{$newChain}++;
            }
            else{print "$line\nSourceCode: " . __LINE__ ."\nFile ($Lfile):$.\n";die}
         }
         else{print "$line\nSourceCode: " . __LINE__ ."\nFile ($Lfile):$.\n";die}
         &end(%finished);

      }#EOF

#end of chain in PISA files           
      elsif( $newChain ne $oldChain   ){

            my ($Fredo) = &newChain($flag,\%finished,\$oldChain,\%temp,$Lfile,$line);
            if($Fredo == 1){redo LINE}
# First finish the old chain
         %temp = ();

         $oldChain = $newChain;
         if($line =~ /^ATOM/){$flag = 1}

         $temp{$newChain}{sequence} .= $resName;
         $temp{$newChain}{records}  .= "$line\n";
         $temp{$newChain}{coord}    .= $x . $y . $z;

         next LINE;
      }
      elsif( $newChain eq $oldChain   ){
         if($line =~ /^ATOM/){$flag = 1}
         $temp{$oldChain}{sequence} .= $resName;
         $temp{$oldChain}{records}  .= "$line\n";
         $temp{$newChain}{coord}    .= $x . $y . $z;

         next LINE;

      }
      else{print "$line\nSourceCode: " . __LINE__ ."\nFile ($Lfile):$.\n";die}
   }    
   sub newChain(){
      if(scalar @_ != 6){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 6 needed\n";die;}
      my ($flag,$HRfinished,$oldChain,$HRtemp,$Lfile,$line) = @_;
      if(   $flag == 0){ 
         if(exists $HRfinished->{onlyHet}->{$$oldChain}){
            if($$oldChain =~ /^\w$/){
               $$oldChain .= '.1';
               return(1);
            }
            elsif(exists $HRfinished->{normal}->{$$oldChain}){
               my @tmp = split(/\./,$$oldChain);
               $tmp[1]++;
               $$oldChain = $tmp[0] . '.' . $tmp[1];
               return(1);
            }
         }
         $HRfinished->{onlyHet}->{$$oldChain} = $HRtemp->{$$oldChain}; 
      }
      elsif($flag == 1){ 
         if(exists $HRfinished->{normal}->{$$oldChain}){
            if($$oldChain =~ /^\w$/){
               $$oldChain .= '.1';
               return(1);
            }
            elsif(exists $HRfinished->{normal}->{$$oldChain}){
               my @tmp = split(/\./,$$oldChain);
               $tmp[1]++;
               $$oldChain = $tmp[0] . '.' . $tmp[1];
               return(1);
            }
         }

         my $originalOldChain = substr($$oldChain,0,1);
         $HRfinished->{normal}->{$$oldChain}  = $HRtemp->{$originalOldChain}; 
         $HRfinished->{chains}->{$$oldChain}++;
      }
      else{print "$line\nSourceCode: " . __LINE__ ."\nFile ($Lfile):$.\n";die}
      return(0);
   }
   sub end(){
      my (%finished) = @_;
# We are at the end of the file, create a few more useful hashes
# Helpful for homooloigomers: Which chain IDs are assigned to a sequence              
         foreach my $chain (sort keys %{$finished{normal}}){
            my $seq = $finished{normal}{$chain}{sequence};
            push(@{$finished{seqToChains}{$seq}},$chain);
#print Dumper(\%finished); 
         }
         foreach my $chain (sort keys %{$finished{normal}}){
            my $coord = $finished{normal}{$chain}{coord};
            $finished{xyzToChains}{$coord} = $chain;
            delete($finished{normal}{$chain}{coord});
#print Dumper(\%finished); 
         }
# Assign one chain as reference chain. Again, useful for homooligomers              
         foreach my $seq ( keys %{$finished{seqToChains}}){
            my $refChain = $finished{seqToChains}{$seq}[0];
            $finished{referenceChain}{$seq} = $refChain; 
            $finished{assignedChains}{$refChain} = @{$finished{seqToChains}{$seq}}; 
         }
         return (\%finished);
   }
}
sub compare_arrays {
   my ($class,$first, $second) = @_;
   no warnings;  # silence spurious -w undef complaints
      return 0 unless @$first == @$second;
   for (my $i = 0; $i < @$first; $i++) {
      return 0 if $first->[$i] ne $second->[$i];
   }
   return 1;
}

sub new {
   my $class = shift;
   my $self = {};
   my %param = @_;
   bless($self, $class);
   return $self;
}
1;

=begin comment      
      if($line eq $lines[$#lines]) {
         my $resName = substr($line,17,3); 
         $resName =~ s/^\s+|\s+$//g; #3
         my $newChain   = substr($line,21,1); 

         if($flag == 0){
            $finished{onlyHet}{$oldChain}{sequence} = $resName;
            $finished{onlyHet}{$oldChain}{records}  = $line;
         }
         elsif($flag == 1){
            $finished{normal}{$oldChain}{sequence} = $resName;
            $finished{normal}{$oldChain}{records}  = $line;
         }
         else{print "$line\nSourceCode: " . __LINE__ ."\nFile ($Lfile):$.\n";die}

         foreach my $chain (sort keys %{$finished{normal}}){
            my $seq = $finished{normal}{$chain}{sequence};
            push(@{$finished{seqToChains}{$seq}},$chain);
#              print Dumper(\%finished); 
         }
# Assign one chain as reference chain. Again, useful for homooligomers              
         foreach my $seq ( keys %{$finished{seqToChains}}){
            my $refChain = $finished{seqToChains}{$seq}[0];
            $finished{referenceChain}{$seq} = $refChain; 
            $finished{assignedChains}{$refChain} = @{$finished{seqToChains}{$seq}}; 
         }
         return (\%finished);
      }
=cut
