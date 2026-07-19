!> Minimal TOML-subset parser for TDSE input files.
!> Supports:
!>   - [section] headers
!>   - [[array_of_tables]] headers
!>   - key = value  (string, integer, float)
!>   - # comments and blank lines
!>   - Quoted strings (single or double)
module toml_parser
    use VarPrecision, only: dp
    implicit none
    private

    ! --- public types and procedures ---
    public :: toml_file
    public :: toml_value
    public :: toml_kv
    public :: toml_table
    public :: toml_error

    !> One key-value pair stored as strings; accessors do type conversions.
    type :: toml_kv
        character(:), allocatable :: section   ! table / array-of-tables name
        character(:), allocatable :: key
        character(:), allocatable :: value     ! raw string from file
        integer :: line  = 0                   ! source line number
        integer :: block = 0                   ! which section-header block this belongs to
    end type toml_kv

    !> Union-like container for a single parsed value (used only internally
    !> for numeric getters, not exposed in the public API beyond toml_file
    !> methods.)
    type :: toml_value
        private
        character(:), allocatable :: raw
    end type toml_value

    !> One table (or array-of-tables element) — essentially a section name
    !> and a list of key-value pairs belonging to it.
    type :: toml_table
        character(:), allocatable :: name
        logical :: is_array = .false.   ! .true. for [[ ... ]]
        type(toml_kv), allocatable :: entries(:)
    end type toml_table

    !> Main parser object — holds all parsed key-value pairs organised by
    !> table (section).  Provides type-specific getter routines.
    type :: toml_file
        type(toml_table), allocatable :: tables(:)
        integer :: num_tables = 0
    contains
        procedure :: parse        => toml_parse
        procedure :: finalise     => toml_finalise
        procedure :: get_string   => toml_get_string
        procedure :: get_int      => toml_get_int
        procedure :: get_real     => toml_get_real
        procedure :: count_array  => toml_count_array_entries
        procedure :: get_array_string => toml_get_array_string
        procedure :: get_array_real   => toml_get_array_real
        procedure :: get_array_int    => toml_get_array_int
    end type toml_file

    !> Simple error-reporting type.
    type :: toml_error
        logical :: flag = .false.
        character(:), allocatable :: message
    end type toml_error

    ! --- internal constants ---
    integer, parameter :: MAX_LINE_LEN = 2048
    integer, parameter :: MAX_ENTRIES  = 500

contains

    ! ------------------------------------------------------------------
    !  toml_parse  —  read and parse a TOML file
    ! ------------------------------------------------------------------
    subroutine toml_parse(this, filename, err)
        class(toml_file), intent(inout) :: this
        character(*), intent(in)      :: filename
        type(toml_error), intent(out) :: err

        integer :: unit, ios, line_no, n_tables, n_entries, current_block
        integer :: entries_this_block
        character(MAX_LINE_LEN) :: line
        character(:), allocatable :: current_section
        logical :: current_is_array, in_table
        type(toml_kv), allocatable :: temp_kvs(:)
        integer :: block_is_array(MAX_ENTRIES)
        character(MAX_LINE_LEN) :: key_str, val_str

        ! initialise
        err%flag = .false.
        current_section = ""          ! top-level (unnamed) table
        current_is_array = .false.
        in_table = .false.
        n_tables = 0
        n_entries = 0
        entries_this_block = 0
        current_block = 1
        block_is_array(:) = 0
        allocate(temp_kvs(MAX_ENTRIES))

        open(newunit=unit, file=adjustl(trim(filename)), status='old', &
             action='read', iostat=ios)
        if (ios /= 0) then
            err%flag = .true.
            err%message = "TOML parser: cannot open file '" // trim(filename) // "'"
            return
        end if

        current_block = 0
        line_no = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line_no = line_no + 1

            call strip_comment(line)
            if (len_trim(line) == 0) cycle    ! blank or comment-only line

            ! --- section header ? ---
            if (is_section_header(line)) then
                ! flush previous section
                if (in_table .and. entries_this_block > 0) then
                    n_tables = n_tables + 1
                end if
                call parse_section_header(line, current_section, current_is_array)
                current_block = current_block + 1
                block_is_array(current_block) = 0
                if (current_is_array) block_is_array(current_block) = 1
                in_table = .true.
                entries_this_block = 0
                cycle
            end if

            ! --- key = value ? ---
            if (index(line, '=') > 0) then
                if (.not. in_table) then
                    ! key-value before first [section] — collect as orphan
                    current_section = ""
                    current_is_array = .false.
                    in_table = .true.
                    current_block = 0
                    entries_this_block = 0
                end if
                n_entries = n_entries + 1
                entries_this_block = entries_this_block + 1
                if (n_entries > MAX_ENTRIES) then
                    err%flag = .true.
                    err%message = "TOML parser: too many key-value entries"
                    close(unit)
                    return
                end if
                call parse_key_value(line, key_str, val_str)
                temp_kvs(n_entries)%section = current_section
                temp_kvs(n_entries)%key     = trim(key_str)
                temp_kvs(n_entries)%value   = trim(val_str)
                temp_kvs(n_entries)%line    = line_no
                temp_kvs(n_entries)%block   = current_block
                cycle
            end if

            ! anything else is an error
            write(err%message, '(A,I0,A)') &
                "TOML parser: unrecognised line ", line_no, ": " // trim(line)
            err%flag = .true.
            close(unit)
            return
        end do
        close(unit)

        ! flush last section
        if (in_table .and. entries_this_block > 0) then
            n_tables = n_tables + 1
        end if

        if (n_tables == 0) then
            err%flag = .true.
            err%message = "TOML parser: no valid entries found in '" // trim(filename) // "'"
            return
        end if

        ! --- reorganise flat kv-list into per-table arrays ---
        call build_tables(this, temp_kvs, n_tables, n_entries, &
                          block_is_array, err)
    end subroutine toml_parse

    ! ------------------------------------------------------------------
    !  Build table list from the flat key-value buffer
    ! ------------------------------------------------------------------
    subroutine build_tables(toml, temp_kvs, n_tables, n_entries, &
                            block_is_array, err)
        type(toml_file), intent(out)   :: toml
        type(toml_kv), intent(in)      :: temp_kvs(:)
        integer, intent(in)            :: n_tables, n_entries
        integer, intent(in)            :: block_is_array(:)
        type(toml_error), intent(out)  :: err

        integer :: i, b, entry_start, entry_end, this_count, n_orphans

        err%flag = .false.
        toml%num_tables = n_tables
        allocate(toml%tables(n_tables))

        ! Warn about orphan entries (block 0, before any [section] header)
        n_orphans = 0
        do i = 1, n_entries
            if (temp_kvs(i)%block == 0) n_orphans = n_orphans + 1
        end do
        if (n_orphans > 0) then
            print*, "TOML WARNING: ", n_orphans, &
                " key-value pair(s) before first [section] header — ignored."
            do i = 1, n_entries
                if (temp_kvs(i)%block == 0) then
                    print*, "  line ", temp_kvs(i)%line, ": ", &
                        trim(temp_kvs(i)%key), " = ", trim(temp_kvs(i)%value)
                end if
            end do
        end if

        ! Build one table per block.  Entries are tagged with block IDs
        ! in the parse loop, and blocks are numbered 1..n_tables.
        entry_start = n_orphans + 1
        do b = 1, n_tables
            ! Find the range of entries belonging to this block.
            entry_end = entry_start - 1
            do i = entry_start, n_entries
                if (temp_kvs(i)%block == b) then
                    entry_end = i
                else if (temp_kvs(i)%block > b) then
                    exit
                end if
            end do

            this_count = entry_end - entry_start + 1
            if (this_count < 1) then
                err%flag = .true.
                err%message = "TOML parser: empty table block"
                return
            end if

            toml%tables(b)%name     = temp_kvs(entry_start)%section
            toml%tables(b)%is_array = (block_is_array(b) == 1)
            allocate(toml%tables(b)%entries(this_count))
            do i = 1, this_count
                toml%tables(b)%entries(i) = temp_kvs(entry_start + i - 1)
            end do

            entry_start = entry_end + 1
        end do
    end subroutine build_tables

    ! ------------------------------------------------------------------
    !  toml_finalise — deallocate internal storage
    ! ------------------------------------------------------------------
    subroutine toml_finalise(this)
        class(toml_file), intent(inout) :: this
        integer :: i
        if (allocated(this%tables)) then
            do i = 1, size(this%tables)
                if (allocated(this%tables(i)%entries)) &
                    deallocate(this%tables(i)%entries)
            end do
            deallocate(this%tables)
        end if
    end subroutine toml_finalise

    ! ------------------------------------------------------------------
    !  Getters for scalar key-value pairs (single-valued keys)
    ! ------------------------------------------------------------------
    function toml_get_string(this, section, key, default) result(val)
        class(toml_file), intent(in) :: this
        character(*), intent(in)     :: section, key, default
        character(:), allocatable    :: val
        integer :: i, j
        val = default
        do i = 1, size(this%tables)
            if (this%tables(i)%name == section) then
                do j = 1, size(this%tables(i)%entries)
                    if (this%tables(i)%entries(j)%key == key) then
                        val = unquote(this%tables(i)%entries(j)%value)
                        return
                    end if
                end do
                return   ! section found, key not — keep default
            end if
        end do
    end function toml_get_string

    function toml_get_int(this, section, key, default) result(val)
        class(toml_file), intent(in) :: this
        character(*), intent(in)     :: section, key
        integer, intent(in)          :: default
        integer                      :: val
        character(:), allocatable    :: raw
        integer :: ios
        raw = this%get_string(section, key, "")
        if (len(raw) == 0) then
            val = default
            return
        end if
        read(raw, *, iostat=ios) val
        if (ios /= 0) val = default
    end function toml_get_int

    function toml_get_real(this, section, key, default) result(val)
        class(toml_file), intent(in) :: this
        character(*), intent(in)     :: section, key
        real(dp), intent(in)         :: default
        real(dp)                     :: val
        character(:), allocatable    :: raw
        integer :: ios
        raw = this%get_string(section, key, "")
        if (len(raw) == 0) then
            val = default
            return
        end if
        read(raw, *, iostat=ios) val
        if (ios /= 0) val = default
    end function toml_get_real

    ! ------------------------------------------------------------------
    !  Array-of-tables helpers
    ! ------------------------------------------------------------------
    function toml_count_array_entries(this, section) result(n)
        class(toml_file), intent(in) :: this
        character(*), intent(in)     :: section
        integer :: n, i
        n = 0
        do i = 1, size(this%tables)
            if (this%tables(i)%name == section .and. &
                this%tables(i)%is_array) then
                n = n + 1
            end if
        end do
    end function toml_count_array_entries

    function toml_get_array_string(this, array_section, idx, key, default) result(val)
        class(toml_file), intent(in) :: this
        character(*), intent(in)     :: array_section, key, default
        integer, intent(in)          :: idx       ! 1-based index into array
        character(:), allocatable    :: val
        integer :: i, cnt, j
        val = default
        cnt = 0
        do i = 1, size(this%tables)
            if (this%tables(i)%name == array_section .and. &
                this%tables(i)%is_array) then
                cnt = cnt + 1
                if (cnt == idx) then
                    do j = 1, size(this%tables(i)%entries)
                        if (this%tables(i)%entries(j)%key == key) then
                            val = unquote(this%tables(i)%entries(j)%value)
                            return
                        end if
                    end do
                    return
                end if
            end if
        end do
    end function toml_get_array_string

    function toml_get_array_real(this, array_section, idx, key, default) result(val)
        class(toml_file), intent(in) :: this
        character(*), intent(in)     :: array_section, key
        integer, intent(in)          :: idx
        real(dp), intent(in)         :: default
        real(dp)                     :: val
        character(:), allocatable    :: raw
        integer :: ios
        raw = this%get_array_string(array_section, idx, key, "")
        if (len(raw) == 0) then
            val = default
            return
        end if
        read(raw, *, iostat=ios) val
        if (ios /= 0) val = default
    end function toml_get_array_real

    function toml_get_array_int(this, array_section, idx, key, default) result(val)
        class(toml_file), intent(in) :: this
        character(*), intent(in)     :: array_section, key
        integer, intent(in)          :: idx, default
        integer                      :: val
        character(:), allocatable    :: raw
        integer :: ios
        raw = this%get_array_string(array_section, idx, key, "")
        if (len(raw) == 0) then
            val = default
            return
        end if
        read(raw, *, iostat=ios) val
        if (ios /= 0) val = default
    end function toml_get_array_int

    ! ==================================================================
    !  PRIVATE helper routines
    ! ==================================================================

    !> Remove everything from the first unquoted '#' onward.
    subroutine strip_comment(line)
        character(*), intent(inout) :: line
        integer :: i
        logical :: in_squote, in_dquote
        in_squote = .false.
        in_dquote = .false.
        do i = 1, len(line)
            if (line(i:i) == "'" .and. .not. in_dquote) in_squote = .not. in_squote
            if (line(i:i) == '"' .and. .not. in_squote) in_dquote = .not. in_dquote
            if (line(i:i) == '#' .and. .not. in_squote .and. .not. in_dquote) then
                line(i:) = ' '
                return
            end if
        end do
    end subroutine strip_comment

    !> True if the trimmed line starts with '['.
    logical function is_section_header(line)
        character(*), intent(in) :: line
        character(:), allocatable :: trimmed
        trimmed = trim(adjustl(line))
        is_section_header = (len(trimmed) > 0 .and. trimmed(1:1) == '[')
    end function is_section_header

    !> Extract section name and detect array-of-tables.
    subroutine parse_section_header(line, section, is_array)
        character(*), intent(in)  :: line
        character(:), allocatable, intent(out) :: section
        logical, intent(out)      :: is_array
        character(:), allocatable :: trimmed
        integer :: i1, i2
        trimmed = trim(adjustl(line))
        if (len(trimmed) >= 4 .and. trimmed(1:2) == '[[' .and. &
            trimmed(len(trimmed)-1:len(trimmed)) == ']]') then
            is_array = .true.
            i1 = 3
            i2 = len(trimmed) - 2
        else if (trimmed(1:1) == '[' .and. &
                 trimmed(len(trimmed):len(trimmed)) == ']') then
            is_array = .false.
            i1 = 2
            i2 = len(trimmed) - 1
        else
            is_array = .false.
            section = trim(trimmed)
            return
        end if
        section = trim(adjustl(trimmed(i1:i2)))
    end subroutine parse_section_header

    !> Split "key = value" into two strings.
    subroutine parse_key_value(line, key, value)
        character(*), intent(in)  :: line
        character(*), intent(out) :: key, value
        integer :: ieq
        ieq = index(line, '=')
        key   = adjustl(line(1:ieq-1))
        value = adjustl(line(ieq+1:))
    end subroutine parse_key_value

    !> Strip surrounding quotes (single or double) from a value string.
    function unquote(str) result(res)
        character(*), intent(in)  :: str
        character(:), allocatable :: res
        integer :: n
        character :: c
        res = trim(adjustl(str))
        n = len(res)
        if (n < 2) return
        c = res(1:1)
        if ((c == '"' .or. c == "'") .and. res(n:n) == c) then
            res = res(2:n-1)
        end if
    end function unquote

end module toml_parser