#!/usr/bin/perl -w
##I created the water file by selecting for all hits that were 17/20 similar
##grep -A7 -B7  "17/20" trout_vs_salmon.water > trout_vs_salmon_17outof20.water

##Update:2.09.2019
#Added a minimum similarity function to always pick the more variable regions.
($water_file, $ref_file) = @ARGV;
$infile_name= $water_file;
$infile_name =~ s/\.water//;
$map_file="$infile_name\_water_keys\.txt";

#Read in the file
open($fh, "<", $water_file) || die "Could not open file $water_file/n $!";
open($map, ">", $map_file);
while (<$fh>) {
    chomp;

    #    $seq{$subj} = "";
    if($_ =~ /\#\s+1:\s+(\S+)/){
        $ref=$1;
    }
    if($_ =~ /\#\s+2:\s+(\S+)/){
        $subj = $1;
        if($subj =~ /(\d+$)/){
            $subj = "$ref\_$1";
        }
        #       $subj=$1;
        #       print "Subject: $subj\n";
    }
    if ($_ =~ /\#\s+Identity:\s+(\d+)\/20/){
        $simil{$subj}=$1;
    }

    if ($_ =~ /\#===/){
        next;
    }

#It is important to use $ref instead of \S+ like I did in the first run as it will get the correct string and not the\
 one beneath as it was happening before..
    if ($_ =~ /^$ref\s+\d+\s+(\S+)/){
        $seq{$subj}=$1;
    }

}

$subj_count=0;

foreach $subj(keys %seq){
    #    $subj_count=$seq{$subj}++;
    ++$subj_count;

}
@sim_values=();

foreach $subj(keys %simil){
    push @sim_values, $simil{$subj};
    @sort_sim = sort { $a <=> $b } @sim_values;
    $min_sim=shift @sort_sim;
}
foreach $subj(keys %simil){
    print $map "Id: $subj\nSequence: $seq{$subj}\nSimilarity: $simil{$subj}\/20\n\n";

}
print $map "All score values :@sim_values \n Minimum value: $min_sim \n";
#foreach $subj(keys %simil){
#    if ($simil{$subj}==$min_sim){
#       print "Hits with highest variability: $ref\n";
#    }
#}
close $map;

print "Finished parsing the water file, $subj_count potential sites will now be aligned to the refence sequence.\n";p\
rint " Make sure that you have used the correct reference (should be the first element on the file title..\n";


open($fh2, "<", $ref_file) || die "Could not open file $ref_file/n $!";

while (<$fh2>) {
    chomp;
    if($_=~ /^(>.+)/){
       $id = $1;
       $seqs{$id} = "";
       next;
   }
    $seqs{$id} .= $_;

    for $subj (keys %seq){
        $idx= index $seqs{$id}, $seq{$subj} ;
        if ($idx > -1){
            $idx_pos{$subj}=$idx + 1;
        }
    }
}
foreach $subj (sort {$idx_pos{$a} <=> $idx_pos{$b}} keys %idx_pos) {
    $offset_pos= $idx_pos{$subj} - 5 ;
    $pam_site=substr $seqs{$id}, $offset_pos, 4;
    if ($simil{$subj} == $min_sim){
#       print "Id: $subj\nLocation:$idx_pos{$subj}\nSimilarity: $simil{$subj}\nPAM: $pam_site\n\n";
    }
}
