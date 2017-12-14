C Append text file src to text file dest, writing a comment from name
      subroutine append_f(dest,src,name)
      implicit none

      integer dest
      integer src
      character*(32) name

      character*128 line

      write(dest,'(a)')'C -----------------------------------------'
      write(dest,'(a)')'C appended file: '//name
      write(dest,'(a)')'C -----------------------------------------'
      
      do while(.true.)
        read(src,'(a)',end=1000)line
        write(dest,'(a)')'C '//line
      enddo
      
 1000 continue ! break line

      write(dest,'(a)')'C -----------------------------------------'
      write(dest,'(a)')'C end file: '//name
      write(dest,'(a)')'C -----------------------------------------'

      return
      end
