# encoding: utf-8
require 'matrix.rb'
require 'mathn.rb'


module SimplexMethodHelper

  def solve func, matrix, b, borders, x = nil
    init_logger
    solve_straight func, matrix, b, borders, x
    make_dual func, matrix, b, borders
  end

  def solve_straight func, matrix, b, borders, x = nil
    func, matrix, b,x = *init_input(func, matrix, b, borders, x)
    x, i_base = *first_stage_simplex_method(func, matrix, b, x, borders)
    straight_simplex_method func, matrix, b, x, i_base, borders
  end

  def make_dual func, matrix, b, borders
    b_dual = func
    func_dual = b
    y = vector_string "y", func_dual.size
    w = vector_string "w", b_dual.size
    v = vector_string "v", b_dual.size
    border_low = borders.map{|i| i[0]}
    border_high = borders.map{|i| i[1]}
    matrix_dual = Matrix.rows(matrix).transpose
    add_messages :make_dual, "#{func_dual} * #{y}' + #{border_high} * #{w}' - #{border_low} * #{v.to_a}' --> min", "#{matrix_dual.to_a} * #{y}' + #{w}'' - #{v}' = #{b_dual}'", "w, v >=0"
  end

  #todo
  def solve_dual func, matrix, b, borders, x = nil
    func, matrix, b,x = *init_input(func, matrix, b, borders, x)
    x, i_base = *first_stage_simplex_method(func, matrix, b, x, borders)
    dual_simplex_method  func, matrix, b, x, i_base, borders
  end

  private

  def vector_string sign, size
    (0...size).map{|i| "#{sign}#{i}"}
  end

  def init_logger
    @messages = {} #:straight => [], :make_dual => [], :dual => []
  end

  def add_messages method_name, *messages
    @messages[method_name] = [] unless @messages[method_name]
    messages.each{|message|@messages[method_name] << message}
  end

  #todo
  def dual_simplex_method  func, matrix, b, x, i_base, borders

  end

  def straight_simplex_method func, matrix, b, x, i_base, borders
    add_messages :straight, "ВТОРАЯ СТАДИЯ", "Граничные условия #{borders}"
    n = 0
    success = false
    until success
      add_messages :straight, "#{n} ИТЕРАЦИЯ ВТОРОЙ СТАДИИ", "Текущее х #{x.to_a}", "Базис #{i_base}"
      x, i_base,success = *iterate(x, i_base, matrix, func, borders)
      add_messages :straight, "Вектор х #{x.to_a} с базисом #{i_base} оптимален!" if  success
      n += 1
    end
    [x, i_base]
  end

  def init_input func, matrix, b, borders, x = nil
    func = Vector.elements func
    matrix = Matrix.rows matrix
    b = Vector.elements b
    x = Vector.elements(x.nil? ? borders.map{|i| i[1]} : x)
    [func, matrix, b,  x]
  end

  def first_stage_matrix matrix,w
    iden = matrix_addition w
    matrix_union matrix, iden
  end

  def first_stage_func a, b
    Vector.elements [0]*a+[-1]*b
  end

  def matrix_addition w
    iden = Matrix.identity(w.size).to_a
    iden = Matrix.rows(iden.each_with_index.map do |i,j|
      i[j] = w[j]>=0? 1: -1
      i
    end)
  end

  def matrix_union matrix1, matrix2
    a = matrix2.to_a
    Matrix.rows(matrix1.to_a.each_with_index.map{|i,j|  i + a[j]})
  end

  def vector_union x, w
    Vector.elements(x.to_a + w.to_a)
  end

  def first_stage_simplex_method func, matrix, b, x, borders
    x_first, matrix_first, func_first, borders_first = *first_stage_init_input(func, matrix, b, x, borders)
    pseudo = (x.size...x_first.size).to_a
    i_base = pseudo
    add_messages :straight, "ПЕРВАЯ СТАДИЯ", "Граничные условия #{borders_first}"
    n = 0
    until pseudo_null? x_first, pseudo
      add_messages :straight, "#{n} ИТЕРАЦИЯ ПЕРВОЙ СТАДИИ", "Текущее значение х #{x_first.to_a}", "Базис #{i_base}"
      x_first, i_base,success = *iterate(x_first, i_base, matrix_first, func_first, borders_first)
      add_messages :straight, "Вектор х #{x_first.to_a} с базисом #{i_base} оптимален!" if  success
      n+=1
    end
    [Vector.elements(x_first[0...x.size]), i_base]
  end

  def pseudo_null? x, pseudo
    pseudo.all?{|i| x[i] == 0}
  end

  def iterate x, i_base, matrix, func, borders
    matrix_base = matrix_base i_base, matrix
    func_base = func_base i_base, func
    u = matrix_base.transpose.inverse*func_base
    i_not_base = (0...func.size).to_a-i_base
    deltas = i_not_base.map{|i| [delta(i, func, matrix, u),i]}
    add_messages :straight, "Базисная матрица #{matrix_base.to_a}", "Базисный целевой вектор #{func_base.to_a}", "Вектор u #{u.to_a}", "Небазисные индексы #{i_not_base}", "Значения дельта #{deltas}"
    delta_i0 = choose_i0 deltas, borders, x
    return [x,i_base,true] if delta_i0.nil?
    l = direction delta_i0, i_not_base, i_base, matrix, matrix_base
    add_messages :straight, "Дельта i0 #{delta_i0}", "Вектор направления L #{l}"
                 tetta0_index = step borders, l, delta_i0[1], x, i_base
    x_new = x+tetta0_index[0]*Vector.elements(l)
    if tetta0_index[1] != delta_i0[1]
      i_base_new = (i_base - [tetta0_index[1]] + [delta_i0[1]]).sort
    else
      i_base_new = i_base
    end
    add_messages :straight, "Шаг Тетта0 #{tetta0_index}", "Новый вектор х #{x_new.to_a}", "Новый базис #{i_base_new}"
    [x_new, i_base_new, false]
  end

  def direction delta_i0, i_not_base, i_base, matrix, matrix_base
    l = [0]*(i_not_base+i_base).size
    l[delta_i0[1]] = delta_i0[0] / delta_i0[0].abs
    lb = (-l[delta_i0[1]]*(matrix_base.inverse*Vector.elements(matrix.transpose.to_a[delta_i0[1]]))).each_with_index{|i,j| l[i_base[j]] = i}
    l
  end

  def step borders, l, i0, x, i_base
    fixnum_max = 2**(0.size * 8 -2) -1
    tetta_i0 = [borders[i0][1]-borders[i0][0],i0]
    tetta = l.each_with_index.map do |i, j|
      if i < 0 && i_base.include?(j)
        [(borders[j][0]-x[j])/i, j]
      elsif i > 0 && i_base.include?(j)
        [(borders[j][1]-x[j])/i, j]
      else
        [fixnum_max,j]
      end
    end
    tetta[i0] =  tetta_i0
    add_messages :straight, "Значения тетта #{tetta.select{|i| i[0]!= fixnum_max}}"
    tetta.min{|a,b| a[0] <=> b[0]}
  end

  def choose_i0 deltas, borders, x
    not_crit_opt = deltas.to_a.select{|delta| delta[0]>0 && x[delta[1]]!=borders[delta[1]][1] || delta[0]<0 && x[delta[1]]!=borders[delta[1]][0]}
    add_messages :straight, "Не выполняется критерий оптимальности для  #{not_crit_opt}"
    not_crit_opt.max{|a,b| a[0].abs<=>b[0].abs}
    #not_crit_opt[0]
  end

  def delta i, func, matrix, u
    func[i] - (Matrix.row_vector(matrix.transpose.to_a[i])*u)[0]
  end

  def matrix_base i_base, matrix
    Matrix.rows(matrix.transpose.to_a.each_with_index.select{|i, j| i_base.include? j}.map{|i| i[0]}).transpose
  end

  def func_base i_base, func
    Vector.elements func.to_a.each_with_index.select{|i,j| i_base.include? j}.map{|i| i[0]}
  end

  def first_stage_init_input func, matrix, b, x, borders
    w = b - matrix*x
    x_first = vector_union x, w.map{|i| i.abs}
    matrix_first = first_stage_matrix matrix,w
    func_first = first_stage_func func.to_a.size, w.to_a.size
    borders_first = borders + w.to_a.map{|i| [0, i.abs]}
    [x_first, matrix_first, func_first, borders_first]
  end

end
