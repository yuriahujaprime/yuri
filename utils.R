### Utilities (functions, constants, etc.) used by the other microbiome R scripts.
### Eric Proffitt

### Prime colors.
primegrey = '#7d7d7d'
primeblue = '#006fe8'
primeviolet = '#8888ff'
primeteal = '#35bec1'
primepink = '#ff5faf'
primeorange = '#ff6400'
primeyellow = '#ffdc4b'
primegreen = '#00c491'
primered = '#f42535'

### Additive log-ratio transform.
alr = function(X, tol=1e-10) {
  Y = log((X[,1:(ncol(X) - 1)] + tol) / (X[,ncol(X)] + tol))
  return(Y)
}

### Centered log-ratio transform.
clr = function(X, tol=1e-10) {
  Y = log((X + tol) / apply((X + tol)^(1 / ncol(X)), 1, prod))
  return(Y)
}

### Inverse centered log-ratio transform.
clrinv = function(X) {
  Y = exp(X)
  Y = Y / rowSums(Y)
  return(Y)
}

### Find cell index in gtable object.
find_cell = function(table, row, col, name="core-fg") {
  table_layout = table$layout
  which(table_layout$t==row & table_layout$l==col & table_layout$name==name)
}

### Function for performing cyclic permutations.
cyclicpermute = function(x, c=1, direction='right') {
  c = c %% length(x)
  
  if (c == 0)
    return(x)
  
  if (direction == 'left')
    x = c(x[-(0:c)], x[0:c])
  
  if (direction == 'right') {
    c = length(x) - c
    x = c(x[-(0:c)], x[0:c])
  }
  
  return(x)
}







