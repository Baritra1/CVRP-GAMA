#If you have a partial solution (i.e. you know a working x and y, this finds the corresponding slack variables)

y=[0,23,23,23,23,24]
d=[0,6,19,12,8]
Q=25
k=18
x=[[0,1,1,1,1,0],[0,0,0,0,0,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[0,0,0,0,0,1]]
n=len(d)

def format_int(to_form,length):
    return list(reversed([int(digit) for digit in bin(to_form)[2:].rjust(length,'0')]))
final=[]
for i in x:
    for j in i:
        final.append(j)
for i in y[1:n]:
    for j in format_int(i,Q.bit_length()):
        final.append(j)

#(13) slack bits
l=0
for i in range(n+1):
    l+=x[0][i]
print("S13",format_int(k-l,k.bit_length()))
for i in format_int(k-l,k.bit_length()):
    final.append(i)
# #(15u) slack bits
# for i in range(1,n):
#     print("S15u",format_int(Q-y[i],Q.bit_length()))
#     for j in format_int(Q-y[i],Q.bit_length()):
#         final.append(j)

# #(15l) slack bits
# for i in range(1,n):
#     print("S15l",format_int(y[i]-d[i],Q.bit_length()))
#     for j in format_int(y[i]-d[i],Q.bit_length()):
#         final.append(j)

#y_j-y_i-(d_j+Q)*x_ij+Q=S
arr=[]
for i in range(n):
    for j in range(n):
        if(i==j):
            continue
        
        # print("S",i,j,y[j]-y[i]-(d[j]+Q)*x[i][j]+Q)
        arr+=format_int(y[j]-y[i]-(d[j]+Q)*x[i][j]+Q,Q.bit_length()+1)
print("S14",arr)
for i in arr:
    final.append(i)
print("full solution",final)