#!/usr/bin/perl -w

#For this script I would like to identify regions on the mitochondria which have specific attributes.
#To begin, I would like regions that, before comparing it to another region:
## - Are around 20 bp long
## - Begin with a few repeated nucleotides such as "TTT"
## - Do not include the above set of nucleotides in the 20 bp long sequence

#When compared to the same region in a different species
## - The promoter region will be different
## - There will be about 2-3 bp differences

##Check if there are any duplicates within the array

($crispr_ref, $ref_seq)= @ARGV;

##ARGV should include: crispr_target_fasta, compare_sp_fasta ;

$crispr_ref = $ARGV[0];
$crispr_name= $crispr_ref =~ s/\.fa//gr;
$ref_seq=$ARGV[1];
$ref_name= $ref_seq =~ s/\.fa//gr;
$map_file= "$crispr_name\_map\.txt";
$grna_file="$crispr_name\_keys\.fa";

open($fh, "<", $crispr_ref) || die "Could not open file $crispr_ref $!/n";

open($map, ">", $map_file);
open($grna, ">", $grna_file) ;

while (<$fh>) {
    chomp;
    if($_=~ /^(>.+)/){
        $id = $1;
        $seqs{$id} = "";
        next;
    }
    $seqs{$id} .= $_;
}

foreach $id(keys %seqs){
    @pams=();
    @guides=();
    if ($seqs{$id} =~ /TTT[CGA]/g){
        @pams= split /TTT[CGA]/ , $seqs{$id};
        @pams=grep { $_ ne '' } @pams;
        $numberpams= scalar @pams;
        @guides= filter_seq(@pams);
        if (scalar !@guides == 0){
            $number_guides= scalar @guides;
        }else{
            print "No gRNA candidates were found here.. \n\n";
        }
        $id_count= 0;
        for $guide(@guides){
            next unless $duplicate{$guide}++;
            print "$guide is duplicated.\n";
        }
        for $guide(@guides){
            $guide_id="";
                        $guide_id= "\>$crispr_name\_".++$id_count;
            $guides_seq{$guide_id}= $guide;
            print $grna "$guide_id\n$guides_seq{$guide_id}\n";


        }
        for $guide_id (keys %guides_seq){
            $idx= index $seqs{$id}, $guides_seq{$guide_id} ;
            if ($idx > -1){
                $idx_pos{$guide_id}=$idx + 1;
            }
        }
        foreach $guide_id (sort {$idx_pos{$a} <=> $idx_pos{$b}} keys %idx_pos) {
            print $map "$guide_id\n$guides_seq{$guide_id}\n$idx_pos{$guide_id}\n";
            }

    }else{
        print "No PAM sites found.. \n\n";
    }

}
close $grna;
close $map;


sub filter_seq {
    @seqs= (@_);
    @f_pam=();
    for $seq (@seqs){
        $pam_sh = substr $seq, 0, 20;
        if (length ($pam_sh) == 20){
            push @f_pam, $pam_sh;
        }
    }
    return @f_pam;

}
close $fh;


##---
##
print "Initiating Waterman-Smith analysis with full penalty scores.\n";

system("./watering_sequencis.sh", $ARGV[1], $grna_file);


##----
$water_file= "$ref_name\_vs\_$crispr_name\.water";

print "Watermnan-Smith Analysis is finished, output file is called $water_file.\nBegin search for most variable regi\
on between the two mitochondrial sequences.\n";

##Preparing output files
$wk_file="$crispr_name\_water_keys\.txt";
$min_file="$crispr_name\_minimum\.txt";

open($fh, "<", $water_file) || die "Could not open file $water_file/n $!";
open($wk, ">", $wk_file);
while (<$fh>) {
#I use the next command here otherwise i get a lot of warning messages for uninitialized values
    next unless $. > 16;
    chomp;
    if($_ =~ /\#\s+1:\s+(\S+)/){
        $ref=$1;
    }
    if($_ =~ /\#\s+2:\s+(\S+)/){
        $subj = $1;
        if($subj =~ /(\d+$)/){
            $subj = "$ref\_$1";
        }
    }
    if ($_ =~ /\#\s+Identity:\s+(\d+)\/20/){
        $simil{$subj}=$1;
    }

    if ($_ =~ /\#===/){
        next;
    }

#It is important to use $ref instead of \S+ like I did in the first run as it will get the correct string and not th\
e one beneath as it was happening before..
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
    print $wk "Id: $subj\nSequence: $seq{$subj}\nSimilarity: $simil{$subj}\/20\n\n";

}
print $wk "All score values :@sim_values \nMinimum value: $min_sim \n";
#print "All score values :@sim_values \nMinimum value: $min_sim \n";
close $wk;

#print "Finished parsing the water file, $subj_count potential sites will now be aligned to the refence sequence.\n"\
;
#print "Make sure that you have used the correct reference (should be the first element on the file title..\n";


open($fh2, "<", $ref_seq) || die "Could not open file $ref_seq/n $!";
open($minim, ">", $min_file);

while (<$fh2>) {
    chomp;
    if($_=~ /^(>.+)/){
        $id = $1;
        #Note to self_ I am using the same variable names as at the start. Maybe it is not that wise to do so in the\
 long run, though I do empty the variable below.
        $seqs{$id} = "";
        next;
    }
    $seqs{$id} .= $_;
#    print "$id\n$seqs{$id}\n";
    for $subj (keys %seq){
        $idx= index $seqs{$id}, $seq{$subj} ;
        if ($idx > -1){
            $idx_pos{$subj}=$idx + 1;
        }
    }
}

foreach $subj (sort {$idx_pos{$a} <=> $idx_pos{$b}} keys %idx_pos) {
    $offset_pos= $idx_pos{$subj} - 5 ;
    $pam_site = substr $seqs{$id}, $offset_pos, 4;
    if ($simil{$subj} == $min_sim){
        print $minim "Id: $subj\nLocation:$idx_pos{$subj}\nSimilarity: $simil{$subj}\nPAM: $pam_site\n\n";
    }
}
close $minim;
