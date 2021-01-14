import json

runs=[1754,1759,1760,1761,1763,1765,1766,1767,1768,1769,1770]
for r in runs:
    jsfile="etc/DYalign_%d.json" % r
    f=open(jsfile)
    params=json.loads(f.read())
    f.close()
    params["general"]["noise"]=1
    fo=open('etc/DYalign_%d_n.json' % r, 'w')
    json.dump(params,fo)
    fo.close()
