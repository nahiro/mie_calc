C  FINDPAR input subroutines.

C----- The basic idea

C  FINDPAR is a poor man's RDPAR.  It implements some of the
C central features of RDPAR --- the ability to find a named set
C of data values in an input file, the relaxation of the need for
C data in the input file to appear in the order in which it's
C used, the provision for comments in the input file --- in
C (almost) standard FORTRAN 77.

C  Among the features of RDPAR that are not implemented in FINDPAR:
C all redirection of input (from one parameter name to another, from
C file to terminal, from file to command line); freer-format data
C specification on input; tracing options.

C----- Using FINDPAR

C  To use FINDPAR, you open your input file using ordinary FORTRAN
C procedures.  You will also use ordinary FORTRAN to read data
C values.  The only difference is that you can use the FINDPAR call
C at any point to position the file to a point after a line containing
C a parameter name.  These parameter names help to document what the
C input values are for; they also don't have to appear in the order in
C which the program uses them.  An input file might look like this:

C ! Parameters for July 16 run.

C NLAYERS
C 5

C LAYER PARAMETERS
C 5,1013,296
C .001,.001,.001,.001,.001,.001,.001,.001

C MOLECULE FLAGS
C t,t,t,t,t,t,t,t

C  By convention, ! or a space introduces a comment line.  A
C parameter name appears on a line of its own, followed by the data
C lines to be read using FORTRAN statements.  Here is a sample of
C FORTRAN code to read this parameter file:

C      OPEN (UNIT = IN_UNIT, FILE = INPUT, STATUS = 'OLD')
C      CALL FINDPAR ( IN_UNIT, 'NLAYERS' )
C      READ (IN_UNIT, *) NLAYERS
C      CALL FINDPAR ( IN_UNIT, 'LAYER PARAMETERS' )
C      READ (IN_UNIT, *) Z, PRESS, TEMP
C      READ (IN_UNIT, *) BRO, CLO, HCHO, NO2, O2, O3, OCLO, SO2
C      CALL FINDPAR ( IN_UNIT, 'MOLECULE FLAGS' )
C      READ (IN_UNIT, *) IF_BRO, IF_CLO, IF_HCHO, IF_NO2, IF_O2, IF_O3,
C     $ IF_OCLO, IF_SO2
C      CLOSE (UNIT = IN_UNIT)

C----- Details of parameter file format

C The parameter file consists of a sequence of comments and
C parameter-name/data chunks, in any order (except as noted below).

C Comments are essentially anything that doesn't happen to get matched
C as a parameter name or read as data.  However, by convention comments
C are introduced by an exclamation point, to make absolutely sure that
C they don't get mistaken for parameter names.  It is also conventional
C to ``comment out'' a parameter name that you don't want to match by
C simply indenting it one space.

C Parameter names can theoretically be almost any sequence of
C characters.  However, FINDPAR ignores trailing spaces both in the
C parameter name you give it and in lines it reads from the parameter
C file.  And, by convention, parameter names don't begin with spaces.
C Parameter names may contain embedded spaces.  Parameter name matching
C is case-insensitive but otherwise exact: if you use odd things
C like control characters in a parameter name they have to match
C exactly, too, so it's usually best to avoid control characters and
C tabs in parameter names.

C The data lines that follow a line with a parameter name are ordinary
C FORTRAN input lines.  FINDPAR has no influence on how they're handled,
C and all the ordinary rules about FORTRAN READ statements apply.

C Comments may not appear on the same lines as parameter names or data,
C and they can't appear within parameter name/data chunks.  When FINDPAR
C locates a parameter name it positions the file at the line following
C that name; your program will then use ordinary FORTRAN READ statements
C to read the data, and will probably be confused by comments that appear
C at that point.

C There's a maximum length for input lines and parameter names, set
C by the parameter MAXLINE in the code.

C By convention, a parameter name that ends in ? precedes an
C input line that contains a single logical variable (T or F).

C----- Repeated parameters

C Sometimes you want to read several parameters with the same name:
C for example, parameters for each layer of the atmosphere, where there
C might be a large number of layers.  In this case the order of
C parameters in the file does matter, since you usually just want to
C put the numbers in the file in some natural order rather than giving
C a number in each input line specifying where it appears in the input
C order.

C The normal FINDPAR approach is to organize these parameters as
C follows: begin with some parameter that appears only once in the file
C (the number of layers, for example); then follow it with the instances
C of the repeated parameter and associated data.  The only-once
C parameter can even be a dummy, with no data following; its importance
C is that reading it gets the parameter file positioned at the right
C point.

C The other problem here is knowing when a list of repeated parameters
C is over.  FINDPAR doesn't have any of the tricks for doing this that
C RDPAR has; you must either know exactly how many instances of the
C repeated parameter there are going to be (by reading some once-only
C parameter that tells you), or else you need some special numbers to
C flag the end of the list (zeros, negative numbers, etc.).  You can't
C just keep looking for the same parameter name indefinitely, because
C FINDPAR will just rewind the file for you and loop through its
C contents forever.

C----- Errors and XFINDPAR

C If FINDPAR can't find a parameter, it prints an error message
C and returns.  Your READ statements following the FINDPAR call are
C very likely to run into errors of their own at this point, since
C you'll be positioned at some random point in the file (usually at
C the beginning, in this version, but that isn't true in all cases).

C If you want to do something more intelligent about such errors,
C you can use XFINDPAR, which is a function returning a logical value:
C .true. if the parameter name was found, .false. if not; XFINDPAR
C doesn't display FINDPAR's error message when a parameter name cannot
C be found.

C Both FINDPAR and XFINDPAR can run into ordinary FORTRAN input errors
C as they read the file: no special action is taken on these---the
C system's default action, whatever that is, occurs.

C       9/10/91         John Lavagnino
C
C  Find a parameter name in a parameter file.

        SUBROUTINE FINDPAR ( NUNIT, PARNAME )

C  Input arguments.

        INTEGER         NUNIT
        CHARACTER * (*) PARNAME

C  Local variables.

        LOGICAL         XFINDPAR
        EXTERNAL        XFINDPAR

C  No local variables need to be SAVEd.

C***********************************************************************

        IF ( .NOT. XFINDPAR ( NUNIT, PARNAME ) ) THEN
          WRITE (*, *) 'FINDPAR error'
          WRITE (*, *) '   Unit number ', NUNIT
          WRITE (*, *) '   Parameter name ', PARNAME
        END IF

        RETURN
        END
C
C  Find a parameter name in a paramter file: return .true. if found,
C .false. if not.

        LOGICAL FUNCTION XFINDPAR ( NUNIT, PARNAME )

C  Parameters.

        INTEGER         MAXLINE
        PARAMETER       ( MAXLINE = 132 )

C  Input arguments.

        INTEGER         NUNIT
        CHARACTER * (*) PARNAME

C  Local variables.

        INTEGER                 III
        INTEGER                 PARLEN
        INTEGER                 START
        INTEGER                 END
        LOGICAL                 ENDSEEN
        CHARACTER * (MAXLINE)   LINE
        CHARACTER * (MAXLINE)   NAME

C  No local variables need to be SAVEd.

C***********************************************************************

C  Determine the length of the parameter name.

        DO 10 III = LEN ( PARNAME ), 1, -1
          IF ( PARNAME ( III : III ) .NE. ' ' ) THEN
            PARLEN = III
            GO TO 20
          END IF
10      CONTINUE

C  If we get to here, then name contains nothing but blanks.
C  We just return, claiming success, in such a case.

        GO TO 500

C  If we get here, then there's a non-null parameter name; but it
C  might still be too long, in which case we always return
C  signaling failure, since we couldn't ever find such a name in the
C  file.

20      CONTINUE
        IF ( PARLEN  .GT.  MAXLINE ) GO TO 400

C  Convert the name to lower-case.

        NAME = PARNAME ( 1 : PARLEN )
        CALL LCSTRING ( NAME ( 1 : PARLEN ) )

C  Top of main loop.

        ENDSEEN = .FALSE.

100     CONTINUE

          LINE = ' '
          READ ( UNIT = NUNIT, FMT = '(A)', END = 200 ) LINE
          CALL LCSTRING ( LINE ( 1 : PARLEN ) )
          IF ( LINE ( 1 : PARLEN ) .NE. NAME ( 1 : PARLEN ) )
     $          GO TO 100

          START = PARLEN + 1
          END = MAXLINE
          IF ( START  .GT.  END ) GO TO 500
          IF ( LINE ( START : END ) .EQ. ' ' ) GO TO 500
          GO TO 100

C  End-of-file branch.

200       CONTINUE
          REWIND ( UNIT = NUNIT )
          IF ( ENDSEEN ) THEN
            GO TO 400
          ELSE
            ENDSEEN = .TRUE.
            GO TO 100
          END IF

C  End of loop: failure, no parameter name found.

400     CONTINUE
        XFINDPAR = .FALSE.
        RETURN

C  End of loop: successful location of parameter name.

500     CONTINUE
        XFINDPAR = .TRUE.
        RETURN

        END
C
C  Find a parameter name in a parameter file, with optional added
C prefix on parameter name.

        SUBROUTINE FINDPAR_PREFIX ( NUNIT, PREFIX, PARNAME )

C  Input arguments.  If PREFIX is something other than a bunch of
C blanks, it will be added onto PARNAME, with a blank to separate
C it, for the FINDPAR call.  If that FINDPAR call fails, then
C PARNAME without prefix is tried instead, and only if that call fails
C is there an error message.  This means that unprefixed versions of
C the parameter name act as defaults.
C  If PREFIX is a bunch of blanks, this is identical to FINDPAR.
C  (This is assuming that PREFIX has no trailing blanks.  Probably
C it should be checking for them and trimming them if necessary.)

        INTEGER         NUNIT
        CHARACTER * (*) PREFIX
        CHARACTER * (*) PARNAME

C  Local variables.

        CHARACTER * 132 LINE

        LOGICAL         XFINDPAR
        EXTERNAL        XFINDPAR

C***********************************************************************

C??? shouldn't do second XFINDPAR if PREFIX is null...

        IF ( PREFIX  .EQ.  ' ' ) THEN
          CALL FINDPAR ( NUNIT, PARNAME )
        ELSE
          LINE = PREFIX // ' ' // PARNAME
          IF ( .NOT. XFINDPAR ( NUNIT, LINE ) ) THEN
            IF ( .NOT. XFINDPAR ( NUNIT, PARNAME ) ) THEN
              WRITE (*, *) 'FINDPAR error'
              WRITE (*, *) '   Unit number ', NUNIT
              WRITE (*, *) '   Parameter name ', PARNAME
              WRITE (*, *) '   Prefix ', PREFIX
            END IF
          END IF
        END IF
        RETURN
        END
C
C  Convert a string to lower-case.  (ASCII character set assumed.)

        SUBROUTINE LCSTRING ( STRING )

C  Modified argument.

        CHARACTER * (*)         STRING

C  Local variables.

        INTEGER         CCC
        INTEGER         III
        INTEGER         OFFSET

C***********************************************************************

        OFFSET = ICHAR ( 'a' ) - ICHAR ( 'A' )
        DO 10 III = 1, LEN ( STRING )
          CCC = ICHAR ( STRING ( III : III ) )
          IF ( CCC .GE. ICHAR ( 'A' ) .AND.
     $          CCC .LE. ICHAR ( 'Z' ) ) THEN
            STRING ( III : III ) = CHAR ( CCC + OFFSET )
          END IF
10      CONTINUE
        RETURN
        END
