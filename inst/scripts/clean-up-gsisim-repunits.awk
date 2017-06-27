

BEGIN {print "repunit collection"}

$1=="REPUNIT" {repu = $2; next}

{print repu, $1}

