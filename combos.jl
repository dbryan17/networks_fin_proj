function bbb(yt, n, w, p) 
  return (yt + n * w * p) / (1 + n * w)
end

k = 0
n = 1000
p = .1
yt = .8
w = 1
new_yt = yt
while k < 10000000
  global new_yt
  global k
  new_yt = bbb(new_yt, n, w, p)
  k += 1
end
print(new_yt)