import irpcana
import json
summary={}
a=irpcana.analyse(1786,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1787,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1788,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1789,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1790,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1791,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1792,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1793,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1794,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1795,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1796,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1797,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1798,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1799,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1800,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1801,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1802,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1803,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1804,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1805,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1806,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1807,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
#a=irpcana.analyse(1808,"etc/DYalign_1786.json");res=a.getratio(rebin=5,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1809,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
a=irpcana.analyse(1810,"etc/DYalign_1786.json");res=a.getratio(rebin=1,ncut=15);summary[res["run"]]=res
#print summary

fout=open("SummaryPR1.json","w");
fout.write(json.dumps(summary, indent=4, sort_keys=True))
fout.close()
