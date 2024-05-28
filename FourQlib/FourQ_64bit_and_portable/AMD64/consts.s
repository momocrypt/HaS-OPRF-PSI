	.file	"consts.c"
	.text
	.globl	PRIME1271
	.data
	.align 32
	.type	PRIME1271, @object
	.size	PRIME1271, 32
PRIME1271:
	.quad	-1
	.quad	9223372036854775807
	.quad	-1
	.quad	9223372036854775807
	.globl	TWOx8
	.align 32
	.type	TWOx8, @object
	.size	TWOx8, 32
TWOx8:
	.long	2
	.long	2
	.long	2
	.long	2
	.long	2
	.long	2
	.long	2
	.long	2
	.globl	ONEx8
	.align 32
	.type	ONEx8, @object
	.size	ONEx8, 32
ONEx8:
	.long	1
	.long	1
	.long	1
	.long	1
	.long	1
	.long	1
	.long	1
	.long	1
	.ident	"GCC: (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	 1f - 0f
	.long	 4f - 1f
	.long	 5
0:
	.string	 "GNU"
1:
	.align 8
	.long	 0xc0000002
	.long	 3f - 2f
2:
	.long	 0x3
3:
	.align 8
4:
