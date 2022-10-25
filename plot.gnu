set size square
set pm3d map

do for [i=0:64] {
  splot './data/u'.i.'.dat' u 1:2:3 title i with pm3d
  pause 0.1
}