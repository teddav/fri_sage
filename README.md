# FRI: A Polynomial Commitment Scheme Explained

You can checkout the associated article at: https://blog.electisec.tech/fri

## Run

Use Sage to run the script:

```bash
sage fri.sage
```

Or use the Docker image for simplicity: https://teddav.github.io/sagemath/

```bash
docker run --rm --platform linux/amd64 -v $(pwd):/app -w /app sagemath/sagemath 'sage ./fri.sage'
```
