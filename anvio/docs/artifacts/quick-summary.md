

This program will summarize the same collection across all of your profile databases. However, the program will use only the first profile database in the argument list to learn about what is in the collection, so it is not exactly necessary to have this collection defined for all of the other profile databases (though one could argue that it is a good idea to do this regardless...). The collection name you provide to this program must be a collection that is present in at least the first profile database in the argument list. 

## Common errors

If you get an error that looks like this:
```
Config Error: The database at [PROFILE.db] does not seem to have a table named
              `detection_splits` :/ Here is a list of table names this database knows:
              [...]
```

That means your profile databases are not the correct version. The tables we are accessing in this program were introduced in profile database version 36. So the solution to this error is to migrate your databases to at least that version, using %{anvi-migrate}s. :)
