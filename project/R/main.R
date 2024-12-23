dyn.load("util.so") # TODO: Corregir NAMESPACE, no importa bien util.so

# Llamar a la función definida en util.so
print_message <- function() {
  .Call("hello_gsl")
}
