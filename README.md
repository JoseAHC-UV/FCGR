# FCGR
Hecho para el proyecto de pregrado "Diagnóstico de cáncer mamario por medio de un enfoque híbrido utilizando representación del juego del caos e inteligencia artificial."

Argumentos de entrada: Archivo.fa, k, -g (opción de generar o no las imágenes), -d (resolución de la imagen en dpi).
Salida: Imagen de referencia con la leyenda y ejes. Imagen de FCGR sola. Archivo .txt con las frecuencias normalizadas separadas por comas.

Requisitos: Solo Biopython para el manejo de la secuencias. Se puede descargar aquí:https://biopython.org/wiki/Download

Uso:

Correr usando Python, por ejemplo con bash:

python ./FCGR.py Archivo.fa -k 5 -g

Esto resultará en una imagen FCGR para 5-mers de todos las secuencias en Archivo.fa concatenadas. 
