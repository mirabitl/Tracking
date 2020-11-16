x=[0,-27.648,-27.648,0,27.648,27.648,0]
y=[0,13.824,41.472,55.296,41.472,13.824,0]

for i in range(6):
  if (x[i]!=x[i+1]):
     a=(y[i+1]-y[i])/(x[i+1]-x[i])
     b=y[i]-a*x[i]
     print "Droite",x[i],x[i+1],a,b
  else:
      print "limit",x[i]
