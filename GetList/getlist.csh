#get_file_list.pl -keys "path,filename" -cond "storage=hpss,production=P11id,trgsetupname=AuAu19_production,filetype=daq_reco_MuDst,filename~st_physics,tof=1" -o Run11_19GeV_P11id.list -delim '/' -limit 0

#get_file_list.pl -keys "path,filename" -cond "storage=hpss,production=P10ik,trgsetupname=AuAu39_production,filetype=daq_reco_MuDst,filename~st_physics,tof=1" -o Run10_39GeV_P10ik.list -delim '/' -limit 0

#get_file_list.pl -keys "path,filename" -cond "storage=hpss,production=P10ik,trgsetupname=AuAu62_production,filetype=daq_reco_MuDst,filename~st_physics,tof=1" -o Run10_62GeV_P10ik.list -delim '/' -limit 0

#get_file_list.pl -keys "path,filename" -cond "storage=hpss,production=P10ik,trgsetupname=AuAu200_production,filetype=daq_reco_MuDst,filename~st_physics,tof=1" -o Run10_200GeV_P10ik.list -delim '/' -limit 0

#get_file_list.pl -keys "path,filename" -cond "storage=nfs,production=P17ii,trgsetupname=AuAu54_production_2017,filetype=daq_reco_MuDst,filename~st_physics,tof=1" -o Run17_54GeV_P17ii_mudst_test.list -delim '/' -limit 10
#get_file_list.pl -keys "path,filename" -cond "storage=nfs,production=P17ii,trgsetupname=AuAu54_production_2017,filetype=daq_reco_PicoDst,filename~st_physics,runnumber[]18161005-18162018,tof=1" -o Run17_54GeV_P17ih_pico.list -delim '/' -limit 0

#get_file_list.pl -keys "path,filename", -delim "/" -cond "production=P17ii,trgsetupname=AuAu54_production_2017,filetype=daq_reco_muDst,filename~st_physics,storage=nfs" -o Run17_54GeV_P17ii_mudst_test.list -limit 10

get_file_list.pl -keys "path,filename" -delim "/" -cond "trgsetupname= production_19GeV_2019,production=P23id,library=SL23d,filetype=daq_reco_picodst,storage=nfs"  -o Run19_19p6GeV_2019_pico.list -limit 0

# get_file_list.pl -keys "path,filename" -delim "/" -cond "production=P18ic,trgsetupname=AuAu54_production_2017,filetype=daq_reco_picoDst,filename~st_physics,storage=local,runnumber=18172010"  -o Run17_54GeV_P18c_pico_embedding_onerun.list -limit 0

#get_file_list.pl -keys "runnumber" -delim "/" -cond "production=P18ic,trgsetupname=AuAu54_production_2017,filetype=daq_reco_picoDst,filename~st_physics,storage=nfs" -o Run17_54GeV_P18c_pico_runnumber.list -limit 0

















































