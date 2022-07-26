This is a multi-talented program that seamlessly updates any anvi'o artifact (i.e., anvi'o databases or config files) to its latest version.

You can provide one or more files and the program will migrate all. However, you must choose whether to migrate the databases safely, or quickly.

If you choose to migrate safely, anvi'o will first make a copy of each file and save it as a backup. In case something goes wrong during the migration, it let you know what happened and will restore your original file from its copy. In the case of a successful migration, the old copy will go away gracefully (and quietly).

This is how you migrate safely:

```
anvi-migrate --migrate-safely *.db
```

Of course, we will always suggest that migrating safely is better, because fewer people get angry at us when we do that. In practice though, making those backup copies takes up extra time and it is unlikely that the migration will fail anyway, so if you have a lot of databases to migrate and are okay with a bit of risk, you have the option to migrate quickly instead. In this case, anvi'o will _not_ copy your databases before starting the migration.

```
anvi-migrate --migrate-quickly *.db
```
Please remember that by living life in the fast lane, you forego your safety net. On the rare occasion that the migration does fail, this program will let you know what happened, leave you with a database that has a `.broken` file extension. In this unlikely event, you can always reach out to us.

### Migrating to a specific version

If your database is a few versions behind the highest available version but for whatever reason you don't want to migrate it all the way, you can specify which version to update your database to. Just use the `-t` flag (note: migrating with this parameter only works on ONE database at a time):

```
anvi-migrate --migrate-safely -t 15 CONTIGS.db
```

Then anvi'o will update your database until it is whatever version you specified and stop. Of course, you cannot provide a version number that is higher than the highest available version. Nor can you provide a number that is lower than your database's current version (ie, backwards migration is not possible).

Not sure what your database's current version is? Try %(anvi-db-info)s.

You can always run `anvi-migrate -v` to learn about the versions of artifact types _your_ anvi'o installation can work with.
