### Run the Docker image

Get a help message from the entrypoint.



Useful flags

- `--rm` deletes the container after is stops running. You can use the command `docker container ls --all` to view stopped containers that have not been deleted.
- `-it` allows for an interactive session.
- `--entrypoint "/bin/bash"` overwrites the entrypoint with the bash binary.
