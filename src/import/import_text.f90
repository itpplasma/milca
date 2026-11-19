module milca_import_text
    ! Key-value view of foreign equilibrium inputs: GVEC INI, TOML and YAML
    ! parameter files, SPECTRE TOML, Fortran namelists (VMEC INDATA, VSEQ,
    ! VMEX MIRROR) and DESC text inputs. Keys are lower case without blanks;
    ! a table entry "(m, n)" of table T becomes "t(m,n)", any other entry
    ! "t_key", so GVEC's INI and TOML spellings coincide.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: kv_store_t, read_lines, lower, normalize_key, parse_reals, &
        parse_logical, unquote, read_gvec_ini, read_toml, read_yaml, &
        read_namelist_group, namelist_groups, mode_key

    type :: kv_item_t
        character(len=:), allocatable :: key, value
    end type

    type :: kv_store_t
        type(kv_item_t), allocatable :: items(:)
        integer :: n = 0
    contains
        procedure :: put => kv_put
        procedure :: find => kv_find
        procedure :: has => kv_has
        procedure :: text => kv_text
        procedure :: real => kv_real
        procedure :: int => kv_int
        procedure :: logical => kv_logical
        procedure :: reals => kv_reals
        procedure :: modes => kv_modes
    end type

contains

    subroutine read_lines(path, lines, ok)
        character(*), intent(in) :: path
        character(len=:), allocatable, intent(out) :: lines(:)
        logical, intent(out) :: ok
        integer :: unit, ios, count, longest, length
        character(len=16384) :: buffer
        ok = .false.
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            allocate (character(1) :: lines(0))
            return
        end if
        count = 0; longest = 1
        do
            read (unit, '(a)', iostat=ios, size=length, advance='no') buffer
            if (is_iostat_end(ios)) exit
            count = count + 1
            longest = max(longest, length)
        end do
        rewind (unit)
        allocate (character(longest) :: lines(count))
        do count = 1, size(lines)
            read (unit, '(a)') lines(count)
        end do
        close (unit)
        ok = .true.
    end subroutine

    pure function lower(text) result(out)
        character(*), intent(in) :: text
        character(len(text)) :: out
        integer :: i, code
        out = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
        end do
    end function

    pure function normalize_key(text) result(out)
        ! lower case, no blanks or quotes, ';' as ','
        character(*), intent(in) :: text
        character(len=:), allocatable :: out
        integer :: i
        out = ''
        do i = 1, len_trim(text)
            select case (text(i:i))
            case (' ', achar(9), '"', "'")
            case (';')
                out = out//','
            case default
                out = out//lower(text(i:i))
            end select
        end do
    end function

    pure function mode_key(name, i, j) result(out)
        character(*), intent(in) :: name
        integer, intent(in) :: i, j
        character(len=:), allocatable :: out
        character(32) :: a, b
        write (a, '(i0)') i
        write (b, '(i0)') j
        out = lower(trim(name))//'('//trim(a)//','//trim(b)//')'
    end function

    pure function unquote(text) result(out)
        character(*), intent(in) :: text
        character(len=:), allocatable :: out
        out = trim(adjustl(text))
        ! a list of strings keeps its quotes
        if (len(out) >= 2) then
            if ((out(1:1) == '"' .or. out(1:1) == "'") .and. out(len(out):len(out)) == out(1:1)) then
                if (index(out(2:len(out) - 1), out(1:1)) == 0) out = out(2:len(out) - 1)
            end if
        end if
        out = trim(adjustl(out))
    end function

    subroutine parse_reals(text, values, ok)
        ! Lists in any of "(/ a, b /)", "[a, b]", "a b c" and "3*0.0".
        character(*), intent(in) :: text
        real(dp), allocatable, intent(out) :: values(:)
        logical, intent(out) :: ok
        character(len=len(text)) :: clean
        character(len=:), allocatable :: token
        integer :: i, start, star, repeat, ios
        real(dp) :: value
        clean = text
        do i = 1, len(clean)
            select case (clean(i:i))
            case ('(', ')', '/', '[', ']', ',', achar(9), '&')
                clean(i:i) = ' '
            case ('d', 'D')
                clean(i:i) = 'e'
            end select
        end do
        allocate (values(0))
        ok = .true.
        i = 1
        do while (i <= len_trim(clean))
            if (clean(i:i) == ' ') then
                i = i + 1
                cycle
            end if
            start = i
            do while (i <= len(clean))
                if (clean(i:i) == ' ') exit
                i = i + 1
            end do
            token = clean(start:i - 1)
            repeat = 1
            star = index(token, '*')
            if (star > 0) then
                read (token(:star - 1), *, iostat=ios) repeat
                if (ios /= 0) then
                    ok = .false.
                    return
                end if
                token = token(star + 1:)
            end if
            read (token, *, iostat=ios) value
            if (ios /= 0) then
                ok = .false.
                return
            end if
            values = [values, spread(value, 1, repeat)]
        end do
    end subroutine

    pure function parse_logical(text, fallback) result(value)
        character(*), intent(in) :: text
        logical, intent(in) :: fallback
        logical :: value
        character(len=:), allocatable :: t
        t = lower(unquote(text))
        value = fallback
        if (len(t) == 0) return
        if (t(1:1) == '.') t = t(2:)
        if (len(t) == 0) return
        if (t(1:1) == 't') value = .true.
        if (t(1:1) == 'f') value = .false.
        if (t == '1') value = .true.
        if (t == '0') value = .false.
    end function

    subroutine kv_put(this, key, value)
        class(kv_store_t), intent(inout) :: this
        character(*), intent(in) :: key, value
        type(kv_item_t), allocatable :: grown(:)
        integer :: at
        character(len=:), allocatable :: k
        k = normalize_key(key)
        at = this%find(k)
        if (at > 0) then
            this%items(at)%value = trim(adjustl(value))
            return
        end if
        if (.not. allocated(this%items)) allocate (this%items(64))
        if (this%n == size(this%items)) then
            allocate (grown(2*this%n))
            grown(:this%n) = this%items(:this%n)
            call move_alloc(grown, this%items)
        end if
        this%n = this%n + 1
        this%items(this%n)%key = k
        this%items(this%n)%value = trim(adjustl(value))
    end subroutine

    integer function kv_find(this, key) result(at)
        class(kv_store_t), intent(in) :: this
        character(*), intent(in) :: key
        character(len=:), allocatable :: k
        k = normalize_key(key)
        do at = 1, this%n
            if (this%items(at)%key == k) return
        end do
        at = 0
    end function

    logical function kv_has(this, key)
        class(kv_store_t), intent(in) :: this
        character(*), intent(in) :: key
        kv_has = this%find(key) > 0
    end function

    function kv_text(this, key, fallback) result(value)
        class(kv_store_t), intent(in) :: this
        character(*), intent(in) :: key, fallback
        character(len=:), allocatable :: value
        integer :: at
        at = this%find(key)
        if (at > 0) then
            value = unquote(this%items(at)%value)
        else
            value = fallback
        end if
    end function

    real(dp) function kv_real(this, key, fallback, ok) result(value)
        class(kv_store_t), intent(in) :: this
        character(*), intent(in) :: key
        real(dp), intent(in) :: fallback
        logical, intent(inout), optional :: ok
        real(dp), allocatable :: values(:)
        logical :: parsed
        value = fallback
        if (.not. this%has(key)) return
        call parse_reals(this%text(key, ''), values, parsed)
        if (parsed .and. size(values) >= 1) then
            value = values(1)
        else if (present(ok)) then
            ok = .false.
        end if
    end function

    integer function kv_int(this, key, fallback) result(value)
        class(kv_store_t), intent(in) :: this
        character(*), intent(in) :: key
        integer, intent(in) :: fallback
        value = nint(this%real(key, real(fallback, dp)))
    end function

    logical function kv_logical(this, key, fallback) result(value)
        class(kv_store_t), intent(in) :: this
        character(*), intent(in) :: key
        logical, intent(in) :: fallback
        value = fallback
        if (this%has(key)) value = parse_logical(this%text(key, ''), fallback)
    end function

    subroutine kv_reals(this, key, values, ok)
        class(kv_store_t), intent(in) :: this
        character(*), intent(in) :: key
        real(dp), allocatable, intent(out) :: values(:)
        logical, intent(out) :: ok
        if (.not. this%has(key)) then
            allocate (values(0))
            ok = .false.
            return
        end if
        call parse_reals(this%text(key, ''), values, ok)
    end subroutine

    subroutine kv_modes(this, name, first, second, values)
        ! every entry name(i,j) with its two integer indices
        class(kv_store_t), intent(in) :: this
        character(*), intent(in) :: name
        integer, allocatable, intent(out) :: first(:), second(:)
        real(dp), allocatable, intent(out) :: values(:)
        character(len=:), allocatable :: prefix, inner
        integer :: at, comma, i, j, ios1, ios2
        real(dp), allocatable :: parsed(:)
        logical :: ok
        prefix = normalize_key(name)//'('
        allocate (first(0), second(0), values(0))
        do at = 1, this%n
            associate (key => this%items(at)%key)
                if (len(key) <= len(prefix)) cycle
                if (key(:len(prefix)) /= prefix) cycle
                if (key(len(key):) /= ')') cycle
                inner = key(len(prefix) + 1:len(key) - 1)
                comma = index(inner, ',')
                if (comma == 0) cycle
                read (inner(:comma - 1), *, iostat=ios1) i
                read (inner(comma + 1:), *, iostat=ios2) j
                if (ios1 /= 0 .or. ios2 /= 0) cycle
                call parse_reals(this%items(at)%value, parsed, ok)
                if (.not. ok .or. size(parsed) < 1) cycle
                first = [first, i]
                second = [second, j]
                values = [values, parsed(1)]
            end associate
        end do
    end subroutine

    pure function strip_comment(line, marks) result(out)
        ! text before the first comment mark outside quotes
        character(*), intent(in) :: line, marks
        character(len=:), allocatable :: out
        integer :: i
        character :: quote
        quote = ' '
        do i = 1, len(line)
            if (quote /= ' ') then
                if (line(i:i) == quote) quote = ' '
            else if (line(i:i) == '"' .or. line(i:i) == "'") then
                quote = line(i:i)
            else if (index(marks, line(i:i)) > 0) then
                out = line(:i - 1)
                return
            end if
        end do
        out = line
    end function

    subroutine read_gvec_ini(path, store, ok)
        ! GVEC INI: key = value, '!' and '#' comments, '&' continues a line
        character(*), intent(in) :: path
        type(kv_store_t), intent(out) :: store
        logical, intent(out) :: ok
        character(len=:), allocatable :: lines(:), line, pending
        integer :: i, eq
        call read_lines(path, lines, ok)
        if (.not. ok) return
        pending = ''
        do i = 1, size(lines)
            line = trim(strip_comment(lines(i), '!#'))
            if (len_trim(line) == 0) cycle
            pending = pending//' '//trim(line)
            if (pending(len_trim(pending):len_trim(pending)) == '&') then
                pending = pending(:len_trim(pending) - 1)
                cycle
            end if
            eq = index(pending, '=')
            if (eq > 1) call store%put(pending(:eq - 1), pending(eq + 1:))
            pending = ''
        end do
    end subroutine

    pure integer function bracket_depth(text) result(depth)
        character(*), intent(in) :: text
        integer :: i
        character :: quote
        depth = 0; quote = ' '
        do i = 1, len(text)
            if (quote /= ' ') then
                if (text(i:i) == quote) quote = ' '
            else if (text(i:i) == '"' .or. text(i:i) == "'") then
                quote = text(i:i)
            else if (text(i:i) == '[' .or. text(i:i) == '{') then
                depth = depth + 1
            else if (text(i:i) == ']' .or. text(i:i) == '}') then
                depth = depth - 1
            end if
        end do
    end function

    pure function table_key(table, key) result(out)
        character(*), intent(in) :: table, key
        character(len=:), allocatable :: out, k
        integer :: i
        k = normalize_key(key)
        out = normalize_key(table)
        do i = 1, len(out)
            if (out(i:i) == '.') out(i:i) = '_'
        end do
        if (len(out) == 0) then
            out = k
        else if (len(k) > 0 .and. k(1:1) == '(') then
            out = out//k
        else
            out = out//'_'//k
        end if
    end function

    subroutine read_toml(path, store, ok)
        ! TOML subset: [tables], key = value, multi-line arrays
        character(*), intent(in) :: path
        type(kv_store_t), intent(out) :: store
        logical, intent(out) :: ok
        character(len=:), allocatable :: lines(:), line, table, pending
        integer :: i, eq
        call read_lines(path, lines, ok)
        if (.not. ok) return
        table = ''; pending = ''
        do i = 1, size(lines)
            line = trim(adjustl(strip_comment(lines(i), '#')))
            if (len(pending) > 0) then
                pending = pending//' '//line
                if (bracket_depth(pending) > 0) cycle
                line = pending
                pending = ''
            end if
            if (len(line) == 0) cycle
            if (line(1:1) == '[' .and. index(line, '=') == 0) then
                table = line(2:index(line, ']', back=.true.) - 1)
                cycle
            end if
            if (bracket_depth(line) > 0) then
                pending = line
                cycle
            end if
            eq = index(line, '=')
            if (eq > 1) call store%put(table_key(table, line(:eq - 1)), line(eq + 1:))
        end do
    end subroutine

    subroutine read_yaml(path, store, ok)
        ! YAML subset: top-level scalars and one level of nested maps; lists
        ! of maps (GVEC's stages) are skipped
        character(*), intent(in) :: path
        type(kv_store_t), intent(out) :: store
        logical, intent(out) :: ok
        character(len=:), allocatable :: lines(:), line, table, key, value
        integer :: i, indent, colon
        logical :: skipping
        call read_lines(path, lines, ok)
        if (.not. ok) return
        table = ''; skipping = .false.
        do i = 1, size(lines)
            line = trim(strip_comment(lines(i), '#'))
            if (len_trim(line) == 0) cycle
            indent = verify(line, ' ') - 1
            line = trim(adjustl(line))
            if (line(1:1) == '-') then
                skipping = .true.
                cycle
            end if
            colon = yaml_colon(line)
            if (colon == 0) cycle
            key = line(:colon - 1)
            value = trim(adjustl(line(colon + 1:)))
            if (indent == 0) then
                skipping = .false.
                if (len(value) == 0) then
                    table = key
                else
                    table = ''
                    call store%put(key, value)
                end if
            else if (.not. skipping .and. len(table) > 0 .and. len(value) > 0) then
                call store%put(table_key(table, key), value)
            end if
        end do
    contains
        pure integer function yaml_colon(text) result(at)
            ! the key's colon; "(0, 1): v" keeps the parenthesised key whole
            character(*), intent(in) :: text
            integer :: start
            start = 1
            if (text(1:1) == '(') start = index(text, ')')
            at = index(text(start:), ':')
            if (at > 0) at = at + start - 1
        end function
    end subroutine

    subroutine namelist_groups(path, names, ok)
        ! names of the &group headers of a Fortran namelist file, lower case
        character(*), intent(in) :: path
        character(32), allocatable, intent(out) :: names(:)
        logical, intent(out) :: ok
        character(len=:), allocatable :: lines(:), line
        integer :: i, stop_at
        call read_lines(path, lines, ok)
        allocate (names(0))
        if (.not. ok) return
        do i = 1, size(lines)
            line = trim(adjustl(strip_comment(lines(i), '!')))
            if (len(line) < 2) cycle
            if (line(1:1) /= '&' .and. line(1:1) /= '$') cycle
            ! the separator after the name sits at line(stop_at + 1)
            stop_at = scan(line(2:), ' ,'//achar(9))
            if (stop_at == 0) stop_at = len(line)
            if (lower(line(2:stop_at)) == 'end') cycle
            names = [character(32) :: names, lower(line(2:stop_at))]
        end do
    end subroutine

    subroutine read_namelist_group(path, group, store, ok)
        ! Lenient reader of one namelist group: NAME(i,j) = values, values may
        ! span lines; unknown names are kept, so no variable list is needed.
        ! A name assigned twice keeps its last non-empty value, as a null
        ! value leaves a namelist variable unchanged.
        character(*), intent(in) :: path, group
        type(kv_store_t), intent(out) :: store
        logical, intent(out) :: ok
        character(len=:), allocatable :: lines(:), line, body, name, value
        integer :: i, at, start, name_end, eq
        logical :: inside, found
        call read_lines(path, lines, ok)
        if (.not. ok) return
        inside = .false.; found = .false.; body = ''
        do i = 1, size(lines)
            line = trim(adjustl(strip_comment(lines(i), '!')))
            if (.not. inside) then
                if (len(line) > len_trim(group) .and. (line(1:1) == '&' .or. line(1:1) == '$')) then
                    if (lower(line(2:len_trim(group) + 1)) == lower(trim(group))) then
                        inside = .true.; found = .true.
                        line = line(len_trim(group) + 2:)
                    end if
                end if
                if (.not. inside) cycle
            end if
            at = group_end(line)
            if (at > 0) then
                body = body//' '//line(:at - 1)
                exit
            end if
            body = body//' '//line
        end do
        ok = found
        if (.not. found) return
        ! split the body at every "name[(...)] =" outside quotes
        name = ''; value = ''
        at = 1
        do while (at <= len(body))
            call assignment_at(body, at, start, name_end, eq)
            if (eq > 0) then
                if (len(name) > 0) call keep(name, value)
                name = body(start:name_end)
                value = ''
                at = eq + 1
            else
                value = value//body(at:at)
                at = at + 1
            end if
        end do
        if (len(name) > 0) call keep(name, value)
    contains
        subroutine keep(key, text)
            character(*), intent(in) :: key, text
            character(len=:), allocatable :: t
            t = trim(adjustl(text))
            do while (len(t) > 0)
                if (t(len(t):len(t)) /= ',') exit
                t = trim(t(:len(t) - 1))
            end do
            if (len(t) == 0 .and. store%has(key)) return
            call store%put(key, t)
        end subroutine
    end subroutine

    pure integer function group_end(line) result(at)
        ! position of the closing '/' or &END outside quotes, 0 if none
        character(*), intent(in) :: line
        integer :: i
        character :: quote
        quote = ' '
        do i = 1, len(line)
            if (quote /= ' ') then
                if (line(i:i) == quote) quote = ' '
            else if (line(i:i) == '"' .or. line(i:i) == "'") then
                quote = line(i:i)
            else if (line(i:i) == '/') then
                at = i
                return
            else if (line(i:i) == '&' .or. line(i:i) == '$') then
                if (i + 2 <= len(line)) then
                    if (lower(line(i + 1:i + 3)) == 'end') then
                        at = i
                        return
                    end if
                end if
            end if
        end do
        at = 0
    end function

    pure subroutine assignment_at(body, at, start, name_end, eq)
        ! an assignment "name[(...)] =" begins at body(at:) after a separator
        character(*), intent(in) :: body
        integer, intent(in) :: at
        integer, intent(out) :: start, name_end, eq
        integer :: i, depth
        eq = 0; start = at; name_end = 0
        if (at > 1) then
            if (index(' ,'//achar(9), body(at - 1:at - 1)) == 0) return
        end if
        if (.not. is_letter(body(at:at))) return
        i = at
        do while (i <= len(body))
            if (.not. (is_letter(body(i:i)) .or. is_digit(body(i:i)) .or. body(i:i) == '_')) exit
            i = i + 1
        end do
        if (i <= len(body)) then
            if (body(i:i) == '(') then
                depth = 0
                do while (i <= len(body))
                    if (body(i:i) == '(') depth = depth + 1
                    if (body(i:i) == ')') depth = depth - 1
                    i = i + 1
                    if (depth == 0) exit
                end do
            end if
        end if
        name_end = i - 1
        do while (i <= len(body))
            if (body(i:i) /= ' ') exit
            i = i + 1
        end do
        if (i > len(body)) return
        if (body(i:i) /= '=') return
        if (i < len(body)) then
            if (body(i + 1:i + 1) == '=') return
        end if
        eq = i
    end subroutine

    pure logical function is_letter(c)
        character, intent(in) :: c
        is_letter = (c >= 'a' .and. c <= 'z') .or. (c >= 'A' .and. c <= 'Z')
    end function

    pure logical function is_digit(c)
        character, intent(in) :: c
        is_digit = c >= '0' .and. c <= '9'
    end function
end module milca_import_text
