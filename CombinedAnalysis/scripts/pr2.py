import irpcana
import json
summary={}
a=irpcana.analyse(1818,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1819,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1820,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1821,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1822,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1823,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1824,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1825,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1826,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1827,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1828,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1829,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1830,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1831,"etc/DYalign_1818.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
fout=open("SummaryPR2.json","w");
fout.write(json.dumps(summary, indent=4, sort_keys=True))
fout.close()
