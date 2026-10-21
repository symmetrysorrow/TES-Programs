! Standalone probe only; this is not PHITS source or a production repair.
program copy_records
  use iso_fortran_env, only: iostat_eor, iostat_end
  implicit none
  character(len=4096) :: chunk
  integer :: ios, n, input_unit, output_unit
  logical :: eof
  character(len=1024) :: infile, outfile
  call get_command_argument(1,infile)
  call get_command_argument(2,outfile)
  open(newunit=input_unit,file=trim(infile),status='old',action='read',form='formatted')
  open(newunit=output_unit,file=trim(outfile),status='replace',action='write',form='formatted')
  eof=.false.
  do while(.not.eof)
    do
      chunk=''
      read(input_unit,'(a)',advance='no',iostat=ios,size=n) chunk
      if(n.gt.0) write(output_unit,'(a)',advance='no') chunk(1:n)
      if(ios.eq.0) cycle
      if(ios.eq.iostat_eor) then
        write(output_unit,'(a)') ''
        exit
      else if(ios.eq.iostat_end) then
        if(n.gt.0) write(output_unit,'(a)') ''
        eof=.true.
        exit
      else
        error stop 'unexpected iostat while copying line'
      endif
    enddo
  enddo
  close(input_unit)
  close(output_unit)
end program copy_records
