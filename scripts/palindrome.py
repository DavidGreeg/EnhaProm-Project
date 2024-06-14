# function to return reverse of a string

def isPalindrome(sequence):
	return sequence == sequence[::-1]

# Example
seq1 = "ATAGAT"
seq2 = "ATAATA"

if isPalindrome(seq2):
	print("yes, seq2: '"+str(seq2)+"' is palindrome")

if not isPalindrome(seq1):
	print("no, seq1: '"+str(seq1)+"' is not palindrome")
