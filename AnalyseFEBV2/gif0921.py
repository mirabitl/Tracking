import hana
dorep=False
res_all=[]
res_all.append(hana.getInfo(74,[7.4],comment="HV Scan  threshold with 10 dac unit mask channels : []",source=0,reprocess=dorep))

res_all.append(hana.getInfo(74,[6,6.2,6.4,6.5,6.6,6.7,6.8,6.9,7.2,7.4,7.5 ],comment="HV Scan  threshold with 10 dac unit mask channels : []",source=0,reprocess=dorep))
res_all.append(hana.getInfo(76,[6,6.2,6.4,6.5,6.6,6.7,6.8,6.9,7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac unit mask channels : []",source=0,reprocess=dorep))
res_all.append(hana.getInfo(77,[7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,15,6,7,8]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(78,[7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,15]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(79, [7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac unit mask channels : [6,7,8]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(88,[7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0]",source=0,reprocess=dorep))

res_all.append(hana.getInfo(81,[7.4],comment="HV Scan  threshold with 12 dac FSM param [0,3,2,2]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(82,[7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac FSM param [1,3,3,2]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(83,[7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac FSM param [1,3,3,3]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(117,[7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac FSM param [1,4,3,3]",source=0,reprocess=dorep))

res_all.append(hana.getInfo(89,[6,6.2,6.4,6.5,6.6,6.7,6.8,6.9,7.2,7.4,7.5,7.6 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,15,6,7,8]",source=220,reprocess=dorep))

res_all.append(hana.getInfo(84,[6,6.2,6.4,6.5,6.6,6.7,6.8,6.9,7.2,7.4,7.5, 7.6 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,15,6,7,8]",source=10,reprocess=dorep))

#,"Priority When Beam ON, Source ON",,,,,,,,
res_all.append(hana.getInfo(119,[6,6.5,6.9,7.2,7.5,7.6 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,15,6,7,8]",source=100,reprocess=dorep))
res_all.append(hana.getInfo(123,[6.9,7.2,7.4,7.6 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0]",source=100,reprocess=dorep))
res_all.append(hana.getInfo(95,[6,6.5,6.9,7.2,7.5,7.6 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,15,6,7,8]",source=46,reprocess=dorep))
res_all.append(hana.getInfo(115,[6, 6.5, 6.9, 7.1, 7.2, 7.5, 7.6 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,15,6,7,8]",source=22,reprocess=dorep))
res_all.append(hana.getInfo(109,[6, 6.5 , 6.9, 7.1, 7.2, 7.3 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,15,6,7,8]",source=10,reprocess=dorep))
res_all.append(hana.getInfo(124,[6 , 6.3, 6.5, 6.7, 6.9, 7.1, 7.3, 7.5, 7.6 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,1,6,7]",source=6.9,reprocess=dorep))
res_all.append(hana.getInfo(129,[6 , 6.3, 6.5, 6.7, 6.9, 7.1, 7.3, 7.5, 7.6 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,1,6,7]",source=10,reprocess=dorep))
res_all.append(hana.getInfo(111,[7.2],comment="HV Scan  threshold with 20 dac unit mask channels : [0,1,8,9]",source=6.9,reprocess=dorep))
res_all.append(hana.getInfo(116,[6, 6.5 , 6.9, 7.1, 7.2, 7.3 ],comment="HV Scan  threshold with 20 dac unit mask channels : [0,1,8,9]",source=6.9,reprocess=dorep))
#,"Priority When Beam ON, Source OFF ",,,,,,,,
res_all.append(hana.getInfo(107,[7.2,7.4,7.5 ],comment="HV Scan  threshold with 20 dac unit mask channels : [0,1,8,9]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(108,[6,6.5,7.0,7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac unit mask channels : [0,1,8,9]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(110,[6,6.5,7.0,7.2,7.4,7.5 ],comment="HV Scan  threshold with 20 dac unit mask channels : [6,7]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(118,[6 ,6.3, 6.5,6.8,7.0,7.1,7.2,7.4,7.5 ],comment="HV Scan  threshold with 12 dac unit mask channels : [6,7]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(121,[6.8,7.0,7.1,7.2,7.4 ],comment="HV Scan  threshold with 12 dac unit mask channels :  [0,1,8,9]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(122,[6.8,7.0,7.1,7.2,7.4 ],comment="HV Scan  threshold with 12 dac unit mask channels :  [0]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(125,[6.8,7.0,7.1,7.2,7.4 ],comment="HV Scan  threshold with 12 dac unit mask channels :  []",source=0,reprocess=dorep))
res_all.append(hana.getInfo(126,[6.8,7.0,7.1,7.2,7.4 ],comment="HV Scan  threshold with 12 dac unit mask channels :  [0,1,6,7]",source=0,reprocess=dorep))
res_all.append(hana.getInfo(127,[7.4],comment="HV Scan threshold with 30 dac unit mask channels : []",source=0,reprocess=dorep))
#res_all.append(hana.getInfo(128,[7.4],comment="HV Scan threshold with 30 dac unit mask channels : [] Power Cycle FEB",source=0,reprocess=dorep))
#print(res_all)
import json
fout=open("summary_gif0921.json","w+")
r=json.dumps(res_all,indent=2)
fout.write(r)
fout.close()
exit()
