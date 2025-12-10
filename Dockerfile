# built on fedora
FROM feiler98/pyomics_fedora

RUN mkdir -p /scratch/tmp/feiler/pyomics_metacells
WORKDIR /scratch/tmp/feiler/pyomics_metacells

COPY . .

# . . means from our computer into our container
RUN pip install --no-cache-dir -r requirements.txt
CMD ["python3", "/home/f/feiler/pyomics_metacells/script.py"]