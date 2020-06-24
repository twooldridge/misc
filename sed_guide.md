Cheat sheet for sed commands, with no particular organization.

### Misc

To view a range of lines, say lines 23 through 25:
```
sed -n '23,25p' $file
```

If you wanted to delete the same lines:
```
sed '23,25d' $file
```

Or say, delete lines starting with the string "dummy"
```
sed '/^dummy/d' $file
```

### Substitions

Starting with the ultra basics, substituting a string, say "dummy", if it appears in a line, with the string "blah"
```
sed 's/dummy/blah/' $file
```

If we wanted to substitute all instances of dummy that appear in the line, we need to tag the end of the command with `g` for global
```
sed 's/dummy/blah/g' $file
```

Or if we wanted to only replace the 2nd occurence, then:
```
sed 's/dummy/blah/2' $file
```
