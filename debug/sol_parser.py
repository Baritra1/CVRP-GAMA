#this parses solutions to be more readable
import sys
import argparse

n=4
k=4
q=24

def parse_and_print(n,k,q,arr):
    a=0
    for i in range(n):
        print(f"x_{i}j",arr[a:a+n+1])
        a+=n+1

    for i in range(1,n):
        print(f"y_{i}",arr[a:a+q.bit_length()])
        a+=q.bit_length()
    print("S_13",arr[a:a+k.bit_length()])
    a+=k.bit_length()
    for i in range(n-1):
        print("S_15u",arr[a:a+q.bit_length()])
        a+=q.bit_length()

    for i in range(n-1):
        print("S_15l",arr[a:a+q.bit_length()])
        a+=q.bit_length()

    for i in range(n*n-n):
        print("S_14",arr[a:a+q.bit_length()+1])
        a+=q.bit_length()+1
    print(arr[a:])
    print(len(arr))

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-n', type=int, help='number of nodes in graph')
    parser.add_argument('-k', type=int, help='number of vehicles in fleet')
    parser.add_argument('-q', type=int, help='max capacity of vehicle')
    parser.add_argument('-arr',type=lambda s: [int(x.strip()) for x in s.split(',') if x.strip()], default=None, help='solution,should be csv enclosed with double quotations')
    
    args = parser.parse_args()
    n=args.n
    k=args.k 
    q=args.q 
    arr=args.arr
    parse_and_print(n,k,q,arr)

# arr = "0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0"
# n = 4
# k = 4
# q = 24
#valid solution for testing