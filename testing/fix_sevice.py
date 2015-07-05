f1 = open('testing_ssyevd.cpp', 'r')
f2 = open('testing_ssyevd1.cpp', 'w')
for line in f1:
	f2.write(line.replace('sevice', 'device'))
f2.close()
f1.close()

f1 = open('testing_sgeev.cpp', 'r')
f2 = open('testing_sgeev1.cpp', 'w')
for line in f1:
	f2.write(line.replace('sevice', 'device'))
f2.close()
f1.close()
