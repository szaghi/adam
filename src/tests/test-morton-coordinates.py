import pymorton as pm

for code in range(20001):
  print(code, pm.deinterleave3(code))
