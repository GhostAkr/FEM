module Quad4Pts

h1(r, s) = 0.25 * (1 + r) * (1 + s)
h2(r, s) = 0.25 * (1 - r) * (1 + s)
h3(r, s) = 0.25 * (1 - r) * (1 - s)
h4(r, s) = 0.25 * (1 + r) * (1 - s)

interFunc = [h1, h2, h3, h4]  # Array of interpolation functions

end  # Quad4Pts