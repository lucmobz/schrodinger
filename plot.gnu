set size square
set pm3d map
do for [i=0:32] {
  splot './data/u'.i.'.dat' u 1:2:3 with pm3d
  pause 0.25
}