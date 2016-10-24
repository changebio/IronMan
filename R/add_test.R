#' Add together two numbers.
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @importClassesFrom DOSE enrichResult
#' @importMethodsFrom DOSE show
#' @importMethodsFrom DOSE summary
#' @importMethodsFrom DOSE plot
#' @importFrom DOSE setReadable
#' @seealso \code{\link{enrichResult-class}}, \code{\link{compareCluster}}
#' @keywords manip
#' @export
#' @examples
#' add(1, 1)
#' add(10, 1)
add_hello <-function(x, y) {
  x +y
}