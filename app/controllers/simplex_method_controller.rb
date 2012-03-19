# encoding: utf-8

class SimplexMethodController < ApplicationController

  include SimplexMethodHelper

  def index
    @title = "Симплекс-метод (все индексы с нуля!!!)"
  end

  def get_messages
    func = params[:func].split(/[\,\+\s]+/).map{|i| i.to_i}
    matrix = params[:matrix].split(/\s*\|\s*/).map{|i| i.split(/[\,\+\s]+/).map{|i| i.to_i}}
    b = params[:b].split(/[\,\+\s]+/).map{|i| i.to_i}
    borders = params[:borders].split(/\s*\|\s*/).map{|i| i.split(/[\,\+\s]+/).map{|i| i.to_i}}
    x = params[:x].split(/[\,\+\s]+/).map{|i| i.to_i} unless params[:x]== ""
    solve func, matrix, b, borders,x
    respond_to do |format|
      format.js
    end
  end

end
