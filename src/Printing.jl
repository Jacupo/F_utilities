

  function Print_matrix(title, matrix)
   figure(title)
   pcolormesh(matrix)
   colorbar()
   ylim(size(matrix,1),0)
  end

  function Print_complex_matrix(title, matrix)
   fig = plt.figure(title,figsize=(11, 5), dpi=80)
   plt.subplots_adjust(wspace=0, hspace=0)
   ax12 = plt.subplot2grid((10,20), (0,0), colspan=10, rowspan=10);
   ax12.set_title("Real part")
   pcolormesh(real.(matrix))
   colorbar()
   ylim(size(matrix,1),0)
   ax2m = plt.subplot2grid((10,20), (0,11), colspan=10, rowspan=10);
   ax2m.set_title("Imaginary part")
   pcolormesh(imag.(matrix))
   colorbar()
   ylim(size(matrix,1),0)
  end
