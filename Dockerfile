# built on fedora
FROM pyomics

RUN mkdir -p pyomics_metacells
WORKDIR ./pyomics_metacells

COPY . .

# . . means from our computer into our container
RUN pip install --no-cache-dir -r requirements.txt
CMD ["python3", "script.py"]