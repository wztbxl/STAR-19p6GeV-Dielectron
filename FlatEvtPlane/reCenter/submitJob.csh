#!/bin/bash
date

dir=$(pwd)
echo $dir

if [ ! -d output_all ]; then
     mkdir -p $dir/output_all
	 ln -s $dir/output_all ./
fi

if [ ! -d script_all ]; then
     mkdir -p $dir/script_all
	 ln -s $dir/script_all ./
fi 

if [ ! -d  log_all ]; then
     mkdir -p $dir/log_all
	 ln -s $dir/log_all ./
fi

rm -rf output_all/*.histo.root
rm -rf script_all/*
rm -rf log_all/*

if [ -e runAll_all.job ]; then
     rm -rf runAll_all.job
fi

echo $PWD

echo "Universe     = vanilla">>runAll_all.job
echo "Notification = never">>runAll_all.job
echo "Requirements = Mips>21000">>runAll_all.job
echo "Initialdir   = $PWD">>runAll_all.job
echo "request_memory = 5000M">>runAll_all.job
echo "GetEnv       = True">>runAll_all.job
echo "+Job_Type    = \"cas\"">>runAll_all.job
echo "">>runAll_all.job
echo "">>runAll_all.job

ifile=0
for FILE in `cat datalist_all`
do
     echo $FILE
     echo "#!/bin/csh">>script_all/all_$ifile.csh
     echo "./reCenter $FILE output_all/$ifile">>script_all/all_$ifile.csh
     chmod 755 script_all/all_$ifile.csh

     echo "Executable       = script_all/all_$ifile.csh">>runAll_all.job
     echo "Output           = log_all/all_$ifile.out">>runAll_all.job
     echo "Error            = log_all/all_$ifile.err">>runAll_all.job
     echo "Log              = log_all/all_$ifile.olog">>runAll_all.job
     echo  "Queue" >>runAll_all.job
     echo  "     " >>runAll_all.job
      
     let "ifile+=1";
done

condor_submit runAll_all.job
