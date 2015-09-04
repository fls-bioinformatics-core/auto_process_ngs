# auto_process auto-completion for bash
#
# This is a 
_auto_process_complete() {
    COMPREPLY=()
    local word=${COMP_WORDS[COMP_CWORD]}
    local completions="$(auto_process.py -h 2>/dev/null | grep -A 100 "^Available commands" | grep -v ^Available | cut -c3- | cut -d' ' -f1 | grep ^${word})"
    COMPREPLY=( $(compgen -W "$completions" -- "$word") )
}

complete -f -F _auto_process_complete auto_process.py
