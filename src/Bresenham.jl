#=
https://github.com/vtavernier/Bresenham.jl/

MIT License

Copyright (c) 2019 Vincent Tavernier

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.=#

#=Changes: 
Included optional args f(x,y,args) in order to pass additional data to f. 
=#
function line(f, x0::Int, y0::Int, x1::Int, y1::Int, args...)
    dx = 2 * (x1 - x0)
    dy = 2 * (y1 - y0)
  
    if dx != 0
      if dy != 0
        if abs(dx) >= abs(dy)
          e = dx / 2
          y = y0
          for x = x0:sign(dx):x1
            f(x, y, args)
            # (a) following code runs even if
            # the loop exits after this iteration
            # but we don't care about that side effect
            e -= sign(dx) * sign(dy) * dy
            if sign(dx) > 0 && e < 0 ||
               sign(dx) < 0 && e >= 0
              y += sign(dy)
              e += dx
            end
          end
        else # abs(dx) < abs(dy) (steep case)
          e = dy / 2
          x = x0
          for y = y0:sign(dy):y1
            f(x, y, args)
            # see (a)
            e -= sign(dx) * sign(dy) * dx
            if sign(dy) > 0 && e < 0 ||
               sign(dy) < 0 && e >= 0
              x += sign(dx)
              e += dy
            end
          end
        end
      else # dy == 0 (horizontal line)
        for x = x0:sign(dx):x1
          f(x, y0, args)
        end
      end
    else # dx == 0 (vertical line)
      for y = y0:sign(dy):y1
        f(x0, y, args)
      end
    end
  end # function