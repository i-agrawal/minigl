#!/usr/bin/perl -w

require strict;

my ($dir, $test_dir, $token) = @ARGV;

open S, "grading-scheme.txt" || die "FAIL: could not open grading scheme.";

my %hashed_tests=();

sub run_test
{
    my ($file,$max_error,$max_time)=@_;
    my $pass_error = 0;
    my $pass_time = 0;

    if(!defined $hashed_tests{$file})
    {
        my $timeout = int(($max_time*1.2*3)/1000)+1;
        if($timeout<2){$timeout=2;}
        `rm -f $dir/file.txt ; cp $test_dir/$file.txt $dir/file.txt`;
        `rm -f $dir/file.png ; cp $test_dir/$file.png $dir/file.png`;
        `cd $dir ; timeout $timeout ./minigl-nogl -i file.txt -s file.png -o $token.txt -n 3 || echo TIMEOUT >> $token.txt`;
        open(my $T, '<', "$dir/$token.txt");
        $_=<$T>;
        my $diff = 1000;
        my $time = $max_time*1.2+1000;
        if(/time: (.*)/)
        {
            $time=$1;
        }
        if(/TIMEOUT/)
        {
            $time="TIMEOUT";
        }
        $_=<$T>;
        if(/diff: (.*)/)
        {
            $diff = $1;
        }
        if(/TIMEOUT/)
        {
            $time="TIMEOUT";
        }
        $hashed_tests{$file}=[$diff,$time];
        close($T);
    }
    my ($diff,$time) = @{$hashed_tests{$file}};
    if($time eq "TIMEOUT")
    {
        print "FAIL: ($file) Test timed out.\n";
    }
    else
    {
        if($diff>$max_error)
        {
            print "FAIL: ($file) Too much error. Actual: $diff  Max: $max_error.\n";
        }
        else
        {
            print "PASS: ($file) diff $diff vs $max_error.\n";
            $pass_error = 1;
        }
        if($time>$max_time)
        {
            print "FAIL: ($file) Too much error. Actual: $time  Max: $max_time.\n";
        }
        else
        {
            print "PASS: ($file) time $time vs $max_time.\n";
            $pass_time = 1;
        }
    }
    return $pass_error && $pass_time;
}

my $total_score=0;
my $max_total_score=0;
while(<S>)
{
    /^\s*\#/ && next;
    /^\s*$/ && next;
    if(/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/)
    {
        my $points=$1;
        my $max_error=$2;
        my $max_time=$3;
        my $file=$4;
        $max_total_score+=$points;

        if(&run_test($file,$max_error,$max_time))
        {
            print "+$points points\n";
            $total_score+=$points;
        }
        else
        {
            print "no points\n";
        }
        
        next;
    }

    print "Unrecognized command: $_";
}

print "FINAL SCORE: $total_score / $max_total_score\n";
