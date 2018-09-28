# Code Sample
# Author: Aaron Sanders
# Simplex Method Linear Programming Solver

import sys
#!/usr/bin/python

text = sys.argv[1] #raw_input("Give input file name: ")

f = open(text, "r")

nums = []

#import data from input file
for word in f.read().split():
	nums.append(int(word))

f.close()

#parse input data for calculations - begin
m = nums[0]
n = nums[1]

#print("m: " + str(m) + "\n")
#print("n: " + str(n) + "\n")



b = []

i = 0
while(i < m):
	b.append(nums[i+2])
	i = i + 1

#print("b: ")
#print(b)
#print("\n")



cT = []
j = 0
while(j < n):
	cT.append(nums[j+m+2])
	j = j + 1

#print("cT: ")
#print(cT)
#print("\n")



startA = m + n + 2;
A = []
z = startA
while(z < len(nums)):
	A.append(nums[z])
	z = z + 1
# end

#print("A: ")
#print(A)
#print("\n")


#slack, identity  matrix
idenM = []
i = 0
j = 0
while(i < m):
	idenM.append(1)
	while(j < m and i < 2):
		idenM.append(0)
		j = j + 1
	i = i + 1
	j = 0
	
#print("slack matrix: ")
#print(idenM)
#print("\n")




sMat = []

i = 0
j = 0
z = 1
k = 0
r = 0
while(i < m):
	while(j < n*z):
		sMat.append(A[j])
		j = j + 1
	while(k < m*z): #adding in slack variables
		sMat.append(idenM[k])
		k = k + 1
	sMat.append(b[i])
	i = i + 1
	z = z + 1
		
i = 0
j = 0
while(i < (n + m + 1)):
	if(i < n):
		sMat.append(-1*cT[i])
	else:
		sMat.append(0)
	i = i + 1

#print("sMat: ")
#print(sMat)
#print("\n")

i = 0
Zs = []
while(i < (n+m)):
	Zs.append(sMat[len(sMat)-(n+m+1)+i])
	i = i + 1 

#print(len(sMat))

#print("Zs: ")
#print(Zs)
#print("\n")

#check to see if all entries non-negative
#
# Possible solutions: bounded-infeasible, unbounded, and bounded-feasible
#
#if everyone is non-negative, the solution is optimal and we are done. 
#  bounded-feasible 
#
#if we come across an undefined value when doing a tableau, 
#  then solution is unbounded
#
#if a slack variable takes on a negative value at completion os Simplex,
#  then the solution is bounded-infeasible

#Bland's Rule - when choosing the pivot column, break any ties
#  by choosing the column with the smallest index

minV = min(Zs)
flag = True
g = 0
index = 0
l = 0
pivotRowStart = 0
pivotRowEnd = 0
Message = ""
while(minV < 0):
	#print("\n")
	#print("\n")
	#print("ITERATION " + str(l+1) + "\n")
	#print("the current matrix: ")
	#print(sMat)
	#print("\n")

	while(g < len(Zs)):
		if(Zs[g] == minV):
			index = g
			break
		g = g + 1


	currCol = []
	i = 0
	while(i < m):
		currCol.append(sMat[index + (i*(n+m+1))])
		i = i + 1

	#print("currCol: ")
	#print(currCol)
	#print("\n")

	ratios = []
	Min = 9999999
	IndexOfMin = 0
	i = 0
	temp = 0
	j = 0
	unbounded = True
	while(i < m):
		if(currCol[i] == 0):
			i = i + 1
			continue

		if(currCol[i] < 0):
			i = i + 1
			continue		

		temp = float(b[i]) / float(currCol[i])
		if(temp < Min and temp > 0):
			Min = temp
		i = i + 1	

	if(Min == 9999999):
		flag = False
		Message = "+inf"
		break;

	while(j < len(currCol)):
		if(currCol[j] != 0):
			if(Min == float(b[j]) / float(currCol[j])): 
				indexR = j
		j = j + 1


	#get pivot row into pivCol - begin
	#print("IndexR: " + str(indexR))
	j = 0
	z = 0
	pivRow = []
	while(j < m):
		if(j == indexR):
			break;
		j = j + 1
	
	k = j * (n + m + 1)

	pivotRowStart = k + z
	#print("j: " + str(j))
	#print("z: " + str(z))
	#print("k: " + str(k))

	#print("pivotRowStart: " + str(pivotRowStart))
	#print((n+m+1)*(j+1))

	if(j == 0):
		while(z < (n+m+1)*(j+1)):
			pivRow.append(sMat[k + z])
			z = z + 1
	elif(j == 1):
		while(z < (n+m+1)*(j)):
			pivRow.append(sMat[k + z])
			z = z + 1
	elif(j > 1):
		while(z < (n+m+1)*(j - 1)):
			pivRow.append(sMat[k + z])
			z = z + 1

	#print("z: " + str(z))
	pivotRowEnd = k + z - 1
	#print("pivotRowEnd: " + str(pivotRowEnd))
	#print(pivRow)

	if(len(pivRow) > 0):
		temp = float(pivRow[index])
		if(pivRow[0] != 1):
			i = 0
			while(i < len(pivRow)):
				pivRow[i] = float(pivRow[i]) / temp
				i = i + 1
	
	#print("pivot range: " + str(pivotRowStart) + " - " + str(pivotRowEnd))
	#print("piv row: ")
	#print(pivRow)
	#print("\n")
	
	#end


	#in this code block can perform Guassian elimination and
	#  create the new sMat 
	q = 0
	p = 0
	w = 0
	i = 0
	nextRow = []
	index2 = index
	pivNum = 0
	nextNum = 0
	temp = 0
	while(i < (m+1)):	
		#print("INDEXR: " + str(indexR))
		#print("INDEX: " + str(index))
		#print("THE NEXTROW BLOCK i: " + str(i))
		pivNum = pivRow[index]
		if(i != indexR):
			k = i * (n + m + 1)
		
			while(q < (n+m+1)): #perform Gaussian elim here
				#print(k + q)
				nextRow.append(sMat[k + q])
				q = q + 1
			#print("INDEX2: " + str(index2))	
			
			#HERE
			if((i-1) == 0):
				#print(nextRow)
				#print(index2 + (n+m+1)*1)
				nextNum = nextRow[index2 + (n+m+1)*(-1)]
			else:
				nextNum = nextRow[index2 + (n+m+1)*(i-1)]
		
			#print("\n")
			#print("NEXT ROW: ")
			#print(nextRow)
			#print("INDEX OF NEXTROW: " + str(index2 + (n+m+1)*(i-1)))
			#print("\n")
	
			#print("index of pivRow: " + str(pivNum))
			#print("index of nextRow: " + str(nextNum))
			#print("nextNum: " + str(nextNum))
			#print("pivNum: " + str(pivNum))
			#print("nextNum / pivNum: " + str(float(nextNum)/float(pivNum)))
			#Guassian elim logic -- return here	
			p = 0
			while(p < (n + m + 1)):
				#print("WTF")
				nextRow[w] = nextRow[w] - (float(nextNum) / float(pivNum))*pivRow[p]
				p = p + 1
				w = w + 1
			p = 0	
	
		
		#print("nextRow (inner loop): ")
		#print(nextRow)
		#print("\n")
		i = i + 1
		p = 0
		q = 0
		#print("current w: " + str(w) + "\n")
	

	if(pivotRowStart == 0):
		u = 0
		while(u < len(pivRow)):
			sMat[u] = pivRow[u]
			u = u + 1
	
		count = 0	
		while(u < (len(nextRow) + n + m + 1)):
			sMat[u] = nextRow[count]
			u = u + 1
			count = count + 1
	else:
		u = 0
		y = 0
		while(u < len(sMat)):
		
			#print("y: " + str(y))
			#print("u: " + str(u) + " len: " + str(len(sMat)))
			if(u >= pivotRowStart and u <= pivotRowEnd):
				f = 0
				#print("Fuck 1")
				pivSpace = pivotRowEnd - pivotRowStart
				#print("pivSpace: " + str(pivSpace))
				#print("fuck")
				while(f < pivSpace):
					#print("f: " + str(f))
					#print("pivRow[f]: " + str(pivRow[f]))
					sMat[u] = pivRow[f]
					u = u + 1
					f = f + 1
				#print("pivRow: ")
				#print(pivRow)
				sMat[u] = pivRow[f]
				u = u + 1
			else:
				#print("Fuck 2")
				#print("nextRow[y]: " + str(nextRow[y]))
				sMat[u] = nextRow[y]
				u = u + 1
				y = y + 1
				

	#print("updated sMat: ")
	#print(sMat)
	#print("\n")
	
	f = 0
	Zs = []
	while(f < (n+m)):
		Zs.append(sMat[len(sMat)-(n+m+1)+f])
		f = f + 1 

	o = 0
	while(o < len(b)):
		b[o] = sMat[(n+m)*(o+1)+o]
		o = o + 1

	#print("FINAL b:")
	#print(b)
	#print("\n")

	#print("Zs (inner loop): ")
	#print(Zs)
	#print("\n")
	
	minV = min(Zs)
	g = 0
	index = 0
	l = l + 1
	i = 0

	#minV = 5
	#print("index of pivot column: " + str(index) + "\n")	

o = open("simplex.out", "w")

i = 0
while(i < len(b)):
	if(b[i] < 0):
		flag = False
		o.write("bounded-infeasible")
		o.write("\n")
		break
	i = i + 1
#cT
if(flag):
	i = 0
	while(i < len(cT)):
		o.write(str(cT[i]))
		o.write(" ")
		i = i + 1
	o.write("\n")
	
	i = 0
	j = i
	count = 1
	z = 0
	startInt = (n + m) * i
	stopInt = (n + m) * (i+1)
	while(i < len(b)):
		k = i
		while(k < len(sMat)):
			if(count > 3):
				break;
			#print("k: " + str(k))
			#print("i: " + str(i))
			#print("count: " + str(count))
			#print("((i+1)*(n+m)) / (n+m) - 1: " + str(((i+1)*(n+m)) / (n + m) - 1))
			if(sMat[k] == 1):
				z = 0
				while(z < len(b)):
					#print("z: " + str(z))
					if(k >= startInt and k <= stopInt):
						o.write(str(b[z]))
						o.write("\n")
						break;
					z = z + 1
					startInt = (n+m)*z	
					stopInt = (n+m)*(z+1)
				#o.write(str(b[(((i+1)*(n + m)) / (n + m)) - 1]))
				#o.write("\n")
				count = 1
				startInt = 0
				stopInt = (n+m)
				break;
			count = count + 1
			k = k + (n + m + i - j + 1)
			#i = i + 1
		i = i + 1
		j = j + 1
		count = 1
		

elif(Message != ""):
	o.write(Message)
	o.write("\n")


o.close()
