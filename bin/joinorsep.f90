!Compile and execute program:
!gfortran filtertest.f90 -o out_filter.x
!./out_filter.x

module variables
	implicit none
	integer :: i,j,io
	integer :: nlines,parnm(2)
	character(30) :: method,FR,FW
	character(200) :: parf,parname(5),parval(3)
	character(:),allocatable :: ids(:)
	integer,allocatable :: snps(:,:)
end module variables

program filter
	use variables
	call GET_COMMAND_ARGUMENT(1,parf)
	open(UNIT=16,FILE=parf)
	do i=1,3
		read(16,*,iostat=io) parname(i),parval(i)
		if (io.ne.0) exit
	end do
	do i=4,5
		read(16,*,iostat=io) parname(i),parnm(i-3)
		if (io.ne.0) exit
	end do
		
	open(UNIT=18,FILE=parval(1))
	open(UNIT=19, FILE= parval(2))
	method=parval(3)
	
	!Count lines in file:
	nlines=0
	do
		read(18,*,iostat=io)
		if (io.ne.0) exit
		nlines=nlines+1
	end do
	rewind(18)
	
	write(*,'(a,a)') 'Input file: ',trim(parval(1))
	write(*,'(a,a)') 'Output file: ',trim(parval(2))
	write(*,'(a,a)') 'Processing operation: ',method
	write(*,'(a,i0)') 'Nlines = ',nlines
	write(*,'(a,i0)') 'ID length = ',parnm(1)
	write(*,'(a,i0)') 'Number of SNPs = ',parnm(2)
	
	!Assign dimensions for ids and snps, as well as the length of their elements:
	allocate( character(len=parnm(1)) :: ids(nlines) )
	allocate( snps(nlines,parnm(2)) )
	
	if (method == 'separate') then
		call separate
	else if  (method == 'join') then
		call join
	end if
end program

subroutine separate
	use variables
	!Read the file:
	write(FR,'("(a",I0,",",I0,"i1)")') parnm(1)+1,parnm(2)
	print *, "Read format: ", FR
	do i=1,nlines
		read(18,FR,iostat=io) ids(i),snps(i,1:parnm(2))
		if (io.ne.0) exit
	end do
	
	!Write the new file (space-separated):
	write(FW,'("(a,a1,i1,",I0,"i2)")') parnm(2)-1
	print *, "Write format: ", FW
	do i=1,nlines
		write(19,FW) trim(ids(i)),' ',snps(i,1),snps(i,2:parnm(2))
	end do
end subroutine separate

subroutine join
	use variables
	!Read the file (space-separated):
	write(FR,'("(a",I0,",",I0,"i2)")') parnm(1),parnm(2)
	print *, "Read format: ", FR
	do i=1,nlines
		read(18,FR,iostat=io) ids(i),snps(i,1:parnm(2))
		if (io.ne.0) exit
	end do
	
	!Write the new file (join SNPs):
	write(FW,'("(a,a1,",I0,"i1)")') parnm(2)
	print *, "Write format: ", FW
	do i=1,nlines
		write(19,FW) ids(i),' ',snps(i,1),snps(i,2:parnm(2))
	end do
end subroutine join
